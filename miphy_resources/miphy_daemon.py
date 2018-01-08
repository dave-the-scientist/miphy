import os, sys, threading
from collections import deque
from random import randint
from flask import Flask, request, render_template, json
from miphy_resources.miphy_instance import MiphyInstance, MiphyValidationError

if sys.version_info >= (3,0): # Python 3.x imports
    from io import StringIO
    from tkinter import Tk as tk_root
    from tkinter.filedialog import asksaveasfilename as saveAs
else: # Python 2.x imports
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO
    from Tkinter import Tk as tk_root
    from tkFileDialog import asksaveasfilename as saveAs

def daemonURL(url):
    return '/daemon' + url

class Daemon(object):
    """Background daemon to serve the miphy requests.

      The local version should be initialized with a server port as an integer,
    and with web_server left as False. The user is expected to call
    new_instance() and process_instance(), open a web browser (if desired), and
    finally call start_server() which runs the flask server, and returns when
    no instances remain.
      The web_server argument should never be set to True on a local version. If
    it is, the user isn't expected to do anything beyond instanciation. The
    program expects that the flask server will be started externally (by uwsgi)
    or similar, and the program isn't meant to ever end.
      This class defines several custom HTTP status codes used to signal errors:
    550 - Specific error validating the user's tree.
    551 - Unknown error validating the user's tree.
    552 - Error parsing options or data sent from web Upload form.
    559 - A request was received with an unrecognized session ID.
    """
    def __init__(self, server_port, web_server=False, instance_timeout_inf=False, verbose=False):
        max_upload_size = 10*1024*1024 # 10 MB
        error_log_lines = 10000
        self.server_port = server_port
        self.web_server = web_server
        self.verbose = verbose
        self.sessions = {} # Holds the miphy instances, with session IDs as keys.
        if not web_server: # Running locally.
            self.sessionID_length = 5 # Length of the unique session ID used.
            self.check_interval = 3 # Repeatedly wait this many seconds between running server tasks.
            self.maintain_wait = 2 # Interval that the page sends a signal to maintain the miphy instance.
            self.allowed_wait = {'after_instance':300, 'page_load':300, 'between_checks':10} # Waits before timing out miphy instances.
            if instance_timeout_inf:
                self.allowed_wait['page_load'] = float('inf')
            self.server_refine_limit = None
            self.server_max_seqs = None
        else: # Live, hosted web server.
            self.sessionID_length = 20
            self.check_interval = 10
            self.maintain_wait = 9
            self.allowed_wait = {'after_instance':120, 'page_load':300, 'between_checks':30}
            self.server_refine_limit = 3000 # Number of sequences above which the server disables the optional refine step
            self.server_max_seqs = 10000 # Max number of sequences for the online version
        # # #  Activity and error logging:
        self.should_quit = threading.Event()
        self.buff_lock = threading.Lock()
        self.log_buffer = StringIO()
        self.error_log = deque([], error_log_lines)
        self.error_occurred = False
        self.signals = [] ## TESTING. i think.
        self.age=0 # TESTING
        # # #  Server setup:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        resources_dir = os.path.join(module_dir, 'resources')
        template_dir = os.path.join(resources_dir, 'templates')
        static_dir = os.path.join(resources_dir, 'static')
        self.server = Flask(__name__, template_folder=template_dir, static_folder=static_dir)
        self.server.config['MAX_CONTENT_LENGTH'] = max_upload_size
        # # #  Server listening routes:
        @self.server.before_first_request
        def setup_tasks():
            if self.web_server: # Setup tasks to start for the web version.
                t = threading.Thread(target=self.web_tasks)
                t.daemon = True; t.start()
            else: # Setup tasks to begin for the local version.
                pass
        # #  Routes used in local version only:
        @self.server.route('/results')
        def render_results():
            return render_template('results.html')
        @self.server.route('/')
        def render_index():
            return render_template('index.html')
        @self.server.route(daemonURL('/save-svg-locally'), methods=['POST'])
        def save_svg():
            svgData = request.form['svgData'].encode('UTF-8')
            root = tk_root()
            filename = saveAs()
            root.destroy()
            if filename:
                if not filename.endswith('.svg'):
                    filename += '.svg'
                with open(filename, 'wb') as f:
                    f.write(svgData)
                ret_msg = 'Svg file saved to %s' % (filename)
            else:
                ret_msg = 'File not saved.'
            return ret_msg
        @self.server.route(daemonURL('/save-csv-locally'), methods=['POST'])
        def save_csv():
            root = tk_root()
            filename = saveAs()
            root.destroy()
            if filename:
                csvStr = request.form['csvStr'].encode('UTF-8')
                if not filename.endswith('.csv'):
                    filename += '.csv'
                with open(filename, 'wb') as f:
                    f.write(csvStr)
                return 'Csv file saved to %s' % (filename)
            else:
                return 'Csv file not saved, as no filename was chosen.'
        # #  Routes used in web version only:
        @self.server.route(daemonURL('/upload-files'), methods=['POST'])
        def upload_files():
            try:
                gene_tree_data = request.files['tree-file'].read()
                info_data = request.files['info-file'].read()
                use_coords = request.form['usecoords']
            except Exception as err:
                return (str(err), 552)
            if use_coords == 'true': use_coords = True
            elif use_coords == 'false': use_coords = False
            else: return ('error parsing usecoords value: "%s"' % use_coords, 552)
            try:
                idnum = self.new_instance(gene_tree_data, info_data, use_coords=use_coords)
            except MiphyValidationError as err:
                return (str(err), 550)
            except Exception as err:
                return (str(err), 551)
            numseqs = self.sessions[idnum].num_sequences
            spc = self.sessions[idnum].species
            action_msg, action_info = '', ''
            if self.server_max_seqs and numseqs > self.server_max_seqs:
                action_msg = 'over seq limit'
                action_info = self.server_max_seqs
            elif use_coords==True and self.server_refine_limit and numseqs > self.server_refine_limit:
                action_msg = 'over refine limit'
                action_info = self.server_refine_limit
            return json.dumps({'idnum':idnum, 'actionmsg':action_msg, 'actioninfo':action_info, 'numseqs':numseqs, 'numspc':len(spc)})
        @self.server.route(daemonURL('/process-data'), methods=['POST'])
        def process_data():
            mi, idnum, msg = self.get_instance()
            if mi == None: return msg
            params = [request.form['ILS'], request.form['dups'], request.form['loss'], request.form['spread']]
            params = tuple(map(float, params))
            mi.processed(params)
            numseqs = mi.num_sequences
            return json.dumps({'numseqs':numseqs, 'numspc':len(mi.species),
                'numclstrs':len(mi.clusters[params])})
        @self.server.route(daemonURL('/get-sessions'), methods=['GET'])
        def get_sessions():
            return json.dumps({'numsessions':len(self.sessions),
                'keys':list(self.sessions.keys()), 'age':self.age})
        # #  Routes used in local and web versions:
        @self.server.route(daemonURL('/get-data'), methods=['POST'])
        def send_data():
            mi, idnum, msg = self.get_instance()
            if mi == None: return msg
            info = {'maintainwait':self.maintain_wait*1000, 'speciestree':mi.species_tree_data,
                'specieslist':mi.species, 'treedata':mi.tree_data, 'initweights':mi.init_weights,
                'sequencenames':mi.sequence_names, 'seqspecies':mi.species_mapping,
                'webversion':self.web_server, 'usecoords':mi.use_coords}
            return json.dumps(info)
        @self.server.route(daemonURL('/cluster-tree'), methods=['POST'])
        def cluster_tree():
            mi, idnum, msg = self.get_instance()
            if mi == None: return msg
            params = (float(request.form['ILS']), float(request.form['dups']),
                    float(request.form['loss']), float(request.form['spread']))
            mi.cluster(params)
            return json.dumps(mi.cluster_list[params])
        @self.server.route(daemonURL('/page-loaded'), methods=['POST'])
        def page_loaded():
            mi, idnum, msg = self.get_instance()
            if mi == None: return msg
            mi.page_loaded()
            return 'page-loaded successful.'
        @self.server.route(daemonURL('/maintain-server'), methods=['POST'])
        def maintain_server():
            mi, idnum, msg = self.get_instance()
            if mi == None: return msg
            mi.maintain()
            return 'maintain-server successful.'
        @self.server.route(daemonURL('/instance-closed'), methods=['POST'])
        def instance_closed():
            mi, idnum, msg = self.get_instance()
            if mi == None: return msg
            del self.sessions[idnum]
            if not self.web_server and len(self.sessions) == 0:
                self.should_quit.set()
            return 'instance-closed successful.'
        @self.server.route('/docs')
        def render_docs():
            return render_template('/docs.html')
        @self.server.route('/contact')
        def render_contact():
            return render_template('/contact.html')
        # # #  TESTING CODE
        @self.server.route('/monitor')
        def render_monitor():
            return render_template('/monitor.html')
        # # #  END OF TESTING.

    def new_instance(self, gene_tree_data, info_data, use_coords=True, coords_file=''):
        if type(info_data) == bytes:
            info_data = info_data.decode()
        if type(gene_tree_data) == bytes:
            gene_tree_data = gene_tree_data.decode()
        idnum = self.generateSessionID()
        self.sessions[idnum] = MiphyInstance(gene_tree_data, info_data, self.allowed_wait, use_coords, coords_file, self.verbose, refine_limit=self.server_refine_limit)
        return idnum
    def process_instance(self, idnum, params):
        self.sessions[idnum].processed(params)
    def start_server(self):
        if self.web_server:
            return False # Only used for local version.
        olderr = sys.stderr
        sys.stderr = self.log_buffer
        t = threading.Thread(target=self.server.run,
            kwargs={'threaded':True, 'port':self.server_port})
        t.daemon = True; t.start()
        try:
            while not self.should_quit.is_set():
                self.should_quit.wait(self.check_interval)
                self.parse_err_logs()
                self.collect_garbage()
            self.parse_err_logs()
            if self.error_occurred:
                print("\nAt least one server call responded with an error. Session log:")
                print(''.join(self.error_log))
        except Exception as error:
            self.parse_err_logs()
            print("\nPython encountered an error. Start of session log:")
            print(''.join(self.error_log))
            print("\nEnd of session log. The error:\n"+str(error))
            # raise
        finally:
            sys.stderr = olderr
    def web_tasks(self):
        if not self.web_server:
            return False # Only used for web version.
        while not self.should_quit.is_set():
            self.should_quit.wait(self.check_interval)
            self.collect_garbage()
    def close(self):
        """Careful with this; the web version should probably never have this
        method actually used."""
        self.should_quit.set()

    # # #  Private methods:
    def get_instance(self, should_fail=False):
        """HTTP status code 559 is used here to indicate a response was requested
        for a session ID that does not exist."""
        if should_fail: return None, 0, ('DEBUG ONLY: Intentional fail.', 588)
        idnum = request.form['session_id']
        if idnum in self.sessions:
            return self.sessions[idnum], idnum, 'session ID is valid.'
        else:
            return None, idnum, ("error, invalid session ID %s." % idnum, 559)
    def generateSessionID(self):
        idnum = ''.join([str(randint(0,9)) for i in range(self.sessionID_length)])
        while idnum in self.sessions:
            idnum = ''.join([str(randint(0,9)) for i in range(self.sessionID_length)])
        return idnum
    def collect_garbage(self):
        self.age += 1 # TESTING
        to_remove = []
        for idnum, mi in self.sessions.items():
            alive = mi.still_alive()
            if not alive:
                to_remove.append(idnum)
        for idnum in to_remove:
            del self.sessions[idnum]
        if not self.web_server: # if personal server with no live instances.
            if len(self.sessions) == 0:
                print('last MIPhy instance closed, shutting down server.')
                self.should_quit.set()
    def parse_err_logs(self):
        with self.buff_lock:
            log_data = self.log_buffer.getvalue()
            self.log_buffer.seek(0)
            self.log_buffer.truncate(0)
        for line in log_data.splitlines(True):
            if '/maintain-server HTTP/1.1" 200' not in line:
                retcode = line.rpartition('-')[0].strip().rpartition('"')[2].strip()
                if retcode not in ('200','304') and '* Running on http://' not in line:
                    self.error_occurred = True
                    print('\nError encountered:\n%s' % line.strip())
                self.error_log.append(line)
