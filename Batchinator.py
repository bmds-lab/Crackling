import tempfile, csv

class Batchinator:
    def __init__(self, batchSize):
        self.workingDir = tempfile.TemporaryDirectory()
        self.currentFile = tempfile.NamedTemporaryFile(mode='w',delete=False,dir=self.workingDir.name)
        self.csvWriter = csv.writer(self.currentFile, delimiter=',', quotechar='"', dialect='unix', quoting=csv.QUOTE_MINIMAL)
        self.batchFiles = []
        self.batchSize = batchSize
        self.entryCount = 0

    def __iter__(self):
        # Close current file
        self.currentFile.close()
        # Record file
        self.batchFiles.append(self.currentFile)
        # yeild the file names
        for file in self.batchFiles:
            yield file.name

    def recordEntry(self, entry):
        # Increase the entry count
        self.entryCount += 1
        # Check if a new file is needed
        if self.entryCount > self.batchSize:
            # Close current file
            self.currentFile.close()
            # Record file
            self.batchFiles.append(self.currentFile)
            # Create new file
            self.currentFile = tempfile.NamedTemporaryFile(mode='w',delete=False,dir=self.workingDir.name)
            # Update csv writer
            self.csvWriter = csv.writer(self.currentFile, delimiter=',', quotechar='"', dialect='unix', quoting=csv.QUOTE_MINIMAL)
            # Reset entry count
            self.entryCount = 1
        # Write entry
        self.csvWriter.writerow(entry)

