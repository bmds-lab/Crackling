import sys

## This class displays, and writes to file, every `print` command    
class Logger(object):
    def __init__(self, outputFile):
        self.terminal = sys.stdout
        self.log = open(outputFile, "w+")

    def __del__(self):
        self.log.close()

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
        
    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        self.terminal.flush()
        self.log.flush()
