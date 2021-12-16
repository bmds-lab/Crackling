import tempfile, csv

class Batchinator:
    def __init__(self, batchSize:int):
        self.workingDir = tempfile.TemporaryDirectory()
        self.currentFile = tempfile.NamedTemporaryFile(mode='w',delete=False,dir=self.workingDir.name)
        self.csvWriter = csv.writer(self.currentFile, delimiter=',', quotechar='"', dialect='unix', quoting=csv.QUOTE_MINIMAL)
        self.batchFiles = []
        self.currentBatch = 0
        self.batchSize = batchSize
        self.entryCount = 0

    def __len__(self):
        return len(self.batchFiles)

    def __iter__(self):
        # Close current file
        self.currentFile.close()
        # Record file
        self.batchFiles.append(self.currentFile)
        # yeild the file names
        for file in self.batchFiles:
            self.currentBatch += 1
            yield file.name

    def recordEntry(self, entry:list):
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

