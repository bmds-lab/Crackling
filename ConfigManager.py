from time import localtime, strftime
from shutil import copyfile
import configparser, os, shutil, sys
import glob

class ConfigManager():
    def __init__(self, filePath, messenger):
        # The configuration
        self._configFilePath = filePath

        # The name of the current configuration
        self._fallbackName = strftime("%Y%m%d%H%M%S", localtime())

        # A list of files that Crackling will process
        self._filesToProcess = []

        # The ConfigParser instance
        self._ConfigParser = configparser.ConfigParser(interpolation=None)

        # A method, passed by reference, to send messages back to the "creator"
        self._sendMsg = messenger

        # Manager setup complete... attempt loading the configuration from file
        self._isConfigured = self._attemptLoadingConfig()
    
        if self._isConfigured:
            self._createListOfFilesToAnalyse()
    
        self._duplicateCracklingCodeAndConfig()
    
        # now that we have initially set things up, we should override __setattr__
        #self.__setattr__ = self._setattr_
    
    def __getitem__(self, arg):
        #if self._isConfigured:
        #    configParserSetItem = self._ConfigParser.__getitem__(arg)
        #    self._isConfigured = self._validateConfig()
        #    if not self._isConfigured:
        #        self._sendMsg('Configuration now invalid!')
        #else:
        #    self._sendMsg('Configuration manager not ready to take values')
        #    exit()
    
        return self._ConfigParser.__getitem__(arg)

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
            
        if success:
            success = self._validateConfig()

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
        
            
            # Ok, everything looks good, lets convert to the new config format
            # We can't use the ConfigParser.read_dict(..) method because the original Dict-formatted
            # config is not formatted correctly to be converted to INI, sigh
            self._ConfigParser.add_section('general')
            
            for firstLayer in CONFIG:
                if isinstance(CONFIG[firstLayer], dict):
                    self._ConfigParser.add_section(firstLayer)
                    for secondLayer in CONFIG[firstLayer]:
                        #if secondLayer == 'delimiter':
                        #    secondLayer = f'"{secondLayer}"'
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
        except Exception as e:
            print(e)
            return False
        #self._ConfigParser['output']['delimiter'] = char(self._ConfigParser['output']['delimiter'])
        return True

    def _validateConfig(self):
        c = self._ConfigParser
        # this method should only be ran once the config has been loaded in
        passed = True
    
        # check the binaries are executable
        for x in [
            c['offtargetscore']['binary'],
            c['bowtie2']['binary'],
            c['rnafold']['binary']
        ]:
            if not shutil.which(x):
                passed = False
                self._sendMsg(f'This binary cannot be executed: {x}')
        
        # check that the 'n' value for the consensus is less than or equal to
        # the number of tools being used
        numToolsInConsesus = self.getNumberToolsInConsensus()
        n = int(c['consensus']['n'])
        
        if n > numToolsInConsesus:
            passed = False
            self._sendMsg(f'The consensus approach is incorrectly set. You have specified {numToolsInConsesus} to be ran but the n-value is {n}. Change n to be <= {numToolsInConsesus}.')
       
        
        c['output']['file'] = os.path.join(c['output']['dir'], f"{self.getConfigName()}-{c['output']['fileName']}")

        if os.path.exists(c['output']['file']):
            passed = False
            self._sendMsg(f"The output file already exists: {c['output']['file']}")
            self._sendMsg(f"To avoid loosing data, please rename your output file.")
            
        return passed

    def _createListOfFilesToAnalyse(self):
        # If the input sequence is a directory:
        if os.path.isdir(self._ConfigParser['input']['exon-sequences']):
            for root, dirs, files in os.walk(self._ConfigParser['input']['exon-sequences']):
                for file in sorted(files, reverse=True):
                    self._filesToProcess.append(os.path.join(self._ConfigParser['input']['exon-sequences'], file))

        # If the input sequence is a file:
        elif os.path.isfile(self._ConfigParser['input']['exon-sequences']):
            self._filesToProcess = [self._ConfigParser['input']['exon-sequences']]
            
        # it's something else
        else:
            self._filesToProcess = glob.glob(self._ConfigParser['input']['exon-sequences'])          
    
    def _duplicateCracklingCodeAndConfig(self):
        pass
    
    def getConfigName(self):
        return self._ConfigParser['general']['name'] or self._fallbackName

    def getNumberToolsInConsensus(self):
        # theres a bug in ConfigParser that makes this messy.
        # it cannot be fixed: https://bugs.python.org/issue10387
        return sum(
            map(lambda x : x == 'True', [
                    self._ConfigParser['consensus']['mm10db'],
                    self._ConfigParser['consensus']['sgrnascorer2'],
                    self._ConfigParser['consensus']['chopchop']
                ]
            )
        )

    def getDatasetSizeBytes(self):
        if self.isConfigured():
            return sum([os.path.getsize(x) for x in self._filesToProcess])

    def isConfigured(self):
        return self._isConfigured

    def getIterFilesToProcess(self):
        fileId = 0
        for file in self._filesToProcess:
            
            c = self._ConfigParser
            
            name = self.getConfigName()
            
            c['rnafold']['input'] = os.path.join(c['output']['dir'], f'{name}-rnafold-input.txt')
            c['rnafold']['output'] = os.path.join(c['output']['dir'], f'{name}-rnafold-output.txt')

            c['offtargetscore']['input'] = os.path.join(c['output']['dir'], f'{name}-{fileId}-offtargetscore-input.txt')
            c['offtargetscore']['output'] = os.path.join(c['output']['dir'], f'{name}-{fileId}-offtargetscore-output.txt')

            c['bowtie2']['input'] = os.path.join(c['output']['dir'], f'{name}-bowtie-input.txt')
            c['bowtie2']['output'] = os.path.join(c['output']['dir'], f'{name}-bowtie-output.txt')

            fileId += 1

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
        
    def getErrLogMethod(self):
        from Logger import Logger
        return Logger(os.path.join(
            self._ConfigParser['output']['dir'], 
            '{}-{}.errlog'.format(
                self._ConfigParser['general']['name'],
                self.getConfigName())
            )
        )

    
if __name__ == '__main__':
    print('Cracking: the configuration manager cannot be called independently.')
    print('Exiting...')
    exit()