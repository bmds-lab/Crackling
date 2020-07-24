from time import localtime, strftime
import configparser, os, shutil, sys

class ConfigManager():
    def __init__(self, filePath, messenger):
        
        # The configuration
        self._configFilePath = filePath

        # Has the configuration manager successfully initialised and is ready to pass to Crackling?
        self._isConfigured = False

        # The name of the current configuration
        self._fallbackName = strftime("%Y%m%d%H%M%S", localtime())

        # A list of files that Crackling will process
        self._filesToProcess = []

        # A list of error/warning messages
        self._messages = []

        # The ConfigParser instance
        self._ConfigParser = configparser.ConfigParser()

        # A method, passed by reference, to send messages back to the "creator"
        self._sendMsg = messenger

        # Setup complete... attempt loading the configuration from file
        self._isConfigured = self._attemptLoadingConfig()
    
        if self._isConfigured:
            self._createListOfFilesToAnalyse()
    
    def __getitem__(self, arg):
        return self._ConfigParser.__getitem__(arg)
    
    def _printer(self, str):
        print(f'ConfigManager: {str}')
    
    def _sendMsg(self, str):
        self._sendMsg(str)
    
    def _attemptLoadingConfig(self):
        filename, fileext = os.path.splitext(self._configFilePath)
        success = False
        
        # Check for v1.0.0 which is a Python dictionary inside a .py file
        # (which will be passed by CLI args without an extension)
        if fileext == '':
            success = self._v1_0_0_to_v1_1_0()
        
        # Check for >v1.0.0, which should now be in the INI format
        if not success:
            # There's some indication that it's INI formatted, but now prove it!
            success = self._read_v1_1_0()

        return success

    def _v1_0_0_to_v1_1_0(self):
        # convert the Python dict config format, which has been deprecated,
        # to the new method using ConfigParser.
        
        # begin by checking if the supplied file is from Version 1.0.0
        try:
            import importlib
            lib = importlib.import_module(self._configFilePath)
            CONFIG = lib.CONFIG
        except:
            self._sendMsg('Yikes!!')
            return False
            
        # we know that the config was designed for v1.0.0 if the user
        # cannot specify which tools the consensus approach uses
        if ({'mm10db', 'sgRNAScorer2', 'CHOPCHOP'} != CONFIG['consensus'].keys()):
            
            # Set the missing config attributes
            CONFIG['consensus']['mm10db'] = True
            CONFIG['consensus']['sgRNAScorer2'] = True
            CONFIG['consensus']['CHOPCHOP'] = True
        
            # make sure that everything else is set
            if not (
                'name'                  in CONFIG                    and
                'consensus'             in CONFIG                    and
                'n'                     in CONFIG['consensus']       and
                'input'                 in CONFIG                    and
                'exon-sequences'        in CONFIG['input']           and
                'offtarget-sites'       in CONFIG['input']           and
                'gff-annotation'        in CONFIG['input']           and
                'bowtie2-index'         in CONFIG['input']           and
                'output'                in CONFIG                    and
                'dir'                   in CONFIG['output']          and
                'fileName'              in CONFIG['output']          and
                'delimiter'             in CONFIG['output']          and
                'offtargetscore'        in CONFIG                    and
                'binary'                in CONFIG['offtargetscore']   and
                'threads'               in CONFIG['offtargetscore']   and
                'score-threshold'       in CONFIG['offtargetscore']   and
                'max-distance'          in CONFIG['offtargetscore']   and
                'sgrnascorer2'          in CONFIG                    and
                'model'                 in CONFIG['sgrnascorer2']    and
                'score-threshold'       in CONFIG['sgrnascorer2']    and
                'bowtie2'               in CONFIG                    and
                'binary'                in CONFIG['bowtie2']         and
                'threads'               in CONFIG['bowtie2']         and
                'rnafold'               in CONFIG                    and
                'binary'                in CONFIG['rnafold']         and
                'threads'               in CONFIG['rnafold']         and
                'low_energy_threshold'  in CONFIG['rnafold']         and
                'high_energy_threshold' in CONFIG['rnafold']         
            ):
                self._sendMsg('Your v1.0.0 configuration is invalid. We suggest updating to the new format, defined as per v1.1.0. See the GitHub repository for a sample configuration file. https://github.com/bmds-lab/Crackling')
                return False
        
            # check the binaries are executable
            passed = True
            for x in [
                CONFIG['offtargetscore']['binary'],
                CONFIG['bowtie2']['binary'],
                CONFIG['rnafold']['binary']
            ]:
                if not shutil.which(x):
                    passed = False
                    self._sendMsg(f'This binary cannot be executed: {x}')
            if not passed:
                return False
        
            # Ok, everything looks good, lets convert to the new config format
            # We can't use the ConfigParser.read_dict(..) method because the original Dict-formatted
            # config is not formatted correctly to be converted to INI, sigh
            self._ConfigParser.add_section('general')
            
            for firstLayer in CONFIG:
                if isinstance(CONFIG[firstLayer], dict):
                    self._ConfigParser.add_section(firstLayer)
                    for secondLayer in CONFIG[firstLayer]:
                        self._ConfigParser.set(firstLayer, secondLayer, str(CONFIG[firstLayer][secondLayer]))
                else:
                    self._ConfigParser.set('general', firstLayer, CONFIG[firstLayer])
         
         
            newConfigFileName = f"{self._configFilePath}.ini"
            self._sendMsg(f'We have transformed your configuration file into the new format. See {newConfigFileName}')
            with open(newConfigFileName, 'w+') as fp:
                self._ConfigParser.write(fp)
        
        return True
            
    def _read_v1_1_0(self):
        try:
            with open(self._configFilePath, 'r') as fp:
                self._ConfigParser.read_file(fp)
        except:
            return False
        return True

    def getConfigName(self):
        if self._isConfigured:
            return self._ConfigParser['general']['name']
        else:
            return self._fallbackName
    
    def _createListOfFilesToAnalyse(self):
        # If the input sequence is a directory:
        if os.path.isdir(self._ConfigParser['input']['exon-sequences']):
            for root, dirs, files in os.walk(self._ConfigParser['input']['exon-sequences']):
                for file in sorted(files, reverse=True):
                    self._filesToProcess.append(file)

        # If the input sequence is a file:
        elif os.path.isfile(self._ConfigParser['input']['exon-sequences']):
            self._filesToProcess = [self._ConfigParser['input']['exon-sequences']]
            
        # it's something else
        else:
            self._filesToProcess = glob.glob(self._ConfigParser['input']['exon-sequences'])          

    def getDatasetSizeBytes(self):
        if self.isConfigured():
            return sum([os.path.getsize(x) for x in self._filesToProcess])
        pass

    def isConfigured(self):
        return self._isConfigured

    def getIterFilesToProcess(self):
        for file in self._filesToProcess:
            
            c = self._ConfigParser
            
            name = self.getConfigName()
            
            c['rnafold']['input'] = os.path.join(c['output']['dir'], f'{name}-rnafold-input.txt')
            c['rnafold']['output'] = os.path.join(c['output']['dir'], f'{name}-rnafold-output.txt')

            c['offtargetscore']['input'] = os.path.join(c['output']['dir'], f'{name}-offtargetscore-input.txt')
            c['offtargetscore']['output'] = os.path.join(c['output']['dir'], f'{name}-offtargetscore-output.txt')

            c['bowtie2']['input'] = os.path.join(c['output']['dir'], f'{name}-bowtie-input.txt')
            c['bowtie2']['output'] = os.path.join(c['output']['dir'], f'{name}-bowtie-output.txt')

            c['output']['file'] = os.path.join(c['output']['dir'], f"{name}-{c['output']['fileName']}")

            yield file

    def getLogMethod(self):
        from Logger import Logger
        return Logger(os.path.join(
            self._ConfigParser['output']['dir'], 
            '{}-{}.log'.format(
                self._ConfigParser['general']['name'],
                self.getConfigName())
            )
        )

if __name__ == '__main__':
    print('Cracking: the configuration manager cannot be called independently.')
    print('Exiting...')
    exit()