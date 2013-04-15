import re
import sys

class ArgParser:
    def __init__(self, cliArgs):
        self.args = cliArgs
        self.params = {}
    
    # use this method to grab REQUIRED command-line parameters
    # flag = the command line flag.  For example: "--model".
    def getArg(self, flag):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            return self.args[i+1]
        else:
            message = "You did not specify the required command-line argument: " + flag 
            print "\n. Error: " + message + "\n"
            #raise AssertionError
            exit(1)
            
    # this is somewhat a hack:
    def setArg(self, flag, content):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            self.args[i] = flag
            self.args[i+1] = content
        else:
            self.args.append(flag)
            self.args.append(content)
            
    # use this method to grab OPTIONAL command-line parameters.
    def getOptionalArg(self, flag):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            return self.args[i+1]
        else:
            return False
    
    def doesContainArg(self, flag):
        return self.args.__contains__(flag)
    
    # use this method to grab OPTIONAL command-line toggles (boolean on/off switches)
    def getOptionalToggle(self, flag):
        if self.args.__contains__(flag):
            return True
        else:
            return False

    # This method will return a list of tokens following 'flag', but not including
    # tokens which start with "--"
    def getList(self, flag, type=str):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            returnList = []
            flagPattern = re.compile("^\-\-.*")
            for j in range( i+1, self.args.__len__() ):
                if re.match(flagPattern, self.args[j] ):
                    return returnList
                else:
                    returnList.append( type(self.args[j]) )
            return returnList
        else:
            print "\n. Hmm, I cannot find the command-line argument ", flag
            exit(1)    

    def getOptionalList(self, flag, type=str):
        #print self.args, flag
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            returnList = []
            flagPattern = re.compile("^\-\-.*")
            for j in range( i+1, self.args.__len__() ):
                if re.match(flagPattern, self.args[j] ):
                    return returnList
                else:
                    returnList.append( type(self.args[j]) )
            return returnList
        else:
            return None
    