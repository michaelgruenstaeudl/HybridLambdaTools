#!/usr/bin/env python2
'''General String Operations'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.06.07.2300"
__status__ = "Working"

#########################
### IMPORT OPERATIONS ###
#########################
# none

###############
### CLASSES ###
###############

class ExciseString:
    ''' Excises string via keyword delimitation.
    Args:
        instring <a>, keyword1 <b>, keyword2 <c>
    Returns:
        outstring
    '''
       
    def __init__(self,a,b,c):
        self.instring = a           # string
        self.keyword1 = b           # string
        self.keyword2 = c           # string
        
    def go(self):
        pos1 = AfterFind(self.instring,self.keyword1).go()
        pos2 = self.instring.find(self.keyword2,pos1)
        return self.instring[pos1:pos2]
        
class ReplaceString:
    ''' Replaces a string via keyword delimitation.
    Args:
        instring <a>, keyword1 <b>, keyword2 <c>, replacement <d>
    Returns:
        original string containing replacement(s)
    Note:
        currently requires that keywords are only found once;
        otherwise use TFL: 
        return self.instring.replace(self.instring[pos1:pos2],self.replacement)
    '''
       
    def __init__(self,a,b,c,d):
        self.instring = a           # string
        self.keyword1 = b           # string
        self.keyword2 = c           # string
        self.replacement = d        # string
        
    def go(self):
        pos1 = AfterFind(self.instring,self.keyword1).go()
        pos2 = self.instring.find(self.keyword2,pos1)
        return self.instring[:pos1]+self.replacement+self.instring[pos2:]
# Legacycode:
#        if self.instring[pos1:pos2]:
#            output = self.instring.replace(self.instring[pos1:pos2],self.replacement)
#        # If string to be replaced were empty:
#        if not self.instring[pos1:pos2]:
#            output = self.instring[:pos1]+self.replacement+self.instring[pos2:]
#        return output

class generatingRelativeValues:
    '''class for generating relative values from a list of absolute values;
       needs inlist <a> as input'''
       
    def __init__(self,a):
        self.inlist = a
        
    def go(self):
        alist = self.inlist[2][self.inlist[2].find("each:")+5:].split(";")  
        summe = sum(float(x.strip()[3:]) for x in alist)
        return "\n".join(str(x)+(str(float(x.strip()[:3])/summe)[:4]) for x in alist)

class ClearSplit:
    ''' Splits string without removing delimiter.
    Args:
        instring <a>, delimiter <b>, rightflag <c>
    Returns:
        outlist
    '''

    def __init__(self,a,b,c):
        self.instring = a
        self.delimiter = b
        self.rightflag = c
        
    def go(self):
        splitList = self.instring.split(self.delimiter)
        if len(splitList) < 2:
            print "*** Error in "+self.__class__.__name__
            print "*** Less than two elements after split."
        outlist = [] 
        if len(splitList) == 2:
            if not self.rightflag:
                outlist = [splitList[0]+self.delimiter, splitList[1]]
            if self.rightflag:            
                outlist = [splitList[0], self.delimiter+splitList[1]]
#        else:
#            outlist = []
#            for a,b in MakePairwise(splitList):
#                if !rightflag:
#                    outlist.append()

        return outlist


class AfterFind:
    '''
    Args:
        instring <a>, keyword <b>
    Returns:
        index number of end of keyword in instring
    '''
    def __init__(self,a,b):
        self.instring = a
        self.keyword = b
                
    def go(self):
        pos1 = self.instring.find(self.keyword)+len(self.keyword)
        return pos1


class RemoveExtension:
    '''class for returning instring without extension (i.e. ".txt" or ".trees";
       needs instring <a> and delimiter <b> as input '''
    def __init__(self,a,b):
        self.instring = a
        self.delimiter = b
                
    def go(self):
        return self.instring[:self.instring.rfind(self.delimiter)]

class is_even:
    '''class for checking if number is even;
       needs number <a> as input '''
    def __init__(self,a):
        self.innumber = a

    def go(self):
        return self.innumber % 2 == 0

class MakePairwise:
    ''' Converts a list into pairs of two.
    Args:
        inlist <a>
    Returns:
        outlist: list of pairs
    '''
    def __init__(self,a):
        self.inlist = a

    def go(self):
        from itertools import tee, izip
        a,b = tee(self.inlist)
        next(b,None)
        return izip(a,b)        

###################
### DEFINITIONS ###
###################

def GenerateRelValues(inlist):
    return generatingRelativeValues(inlist).go()
    
def afind(instring,keyword):
    return AfterFind(instring,keyword).go()

def csplit(instring,delimiter,rightflag=False):
    return ClearSplit(instring,delimiter,rightflag).go()

def exstr(instring,keyword1,keyword2):
    return ExciseString(instring,keyword1,keyword2).go()

def replstr(instring,keyword1,keyword2,replacement):
    return ReplaceString(instring,keyword1,keyword2,replacement).go()

def rmext(a):
    return RemoveExtension(a,".").go()

def iseven(a):
    return is_even(a).go()

