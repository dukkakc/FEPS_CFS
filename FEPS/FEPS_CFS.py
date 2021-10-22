#from lib.phychemProbs import *
#from lib.ProteinDescriptors import *
#from django.conf import settings
import os, sys, getopt
from itertools import product
from itertools import combinations
from Bio import SeqIO
from propy import ProCheck
#import os
import string

#from lib.aaindex import *
#import aaindex

#import sys
#import os
import re
from cogent3.parse.record_finder import DelimitedRecordFinder
#from string import rstrip
#from cogent3.maths.matrix.distance import DistanceMatrix
from cogent3.evolve.fast_distance import DistanceMatrix



class AAIndexParser(object):
    def __init__(self):
        """ Initialize the object. """
        
    def __call__(self, infile):
        """ Parse AAIndex file into dict of AAIndex objects with ID as key

            infile = file to parse as file object or list of lines

            Usage:
                aa1p = AAIndex1Parser()
                aaIndex1Objects = aa1p('data/AAIndex1')

                aa2p = AAIndex2Parser()
                aaIndex2Objects = aa2p('data/AAIndex2')
        """
        
        result = {}

        # Break down the file into records delimited by '//' and then
        # parse each record into AAIndexRecord objects which will be stored
        # in a dict keyed by the records unique ID string
        AAIndexRecordFinder = DelimitedRecordFinder('//', constructor=rstrip)
        # parser is a generator of AAIndexRecords from file
        parser = AAIndexRecordFinder(infile)       

        for r in parser:
            new_record = self._parse_record(r)
            if new_record:
                yield new_record

    def _get_field(self, field_identifier, lines):
        """ Returns the field identified as a one line string
        """
        i = 0
        result = ''
        # Concatenate multi-line data with line_split
        line_split = ' '
        # Run through all lines in the current record
        while (i < len(lines)):
            # Check each line to see if it starts with the field
            # identifier we are looking for
            if (lines[i].startswith(field_identifier)):
                # If we find the line we are looking for, include it in
                # the result, unless it's a Data line.
                # Data entries are multi-line, and the first is information
                # that we are not interested in here.
                if (field_identifier != 'I'):
                    result += lines[i]
                    if field_identifier == 'M': result += 'BRK'
                    # Get rid of the line identifier and leading white space
                    result = result[2:]
                # Move to next line
                i += 1
                # and see if it's a continuation from the above line
                while (i < len(lines) and\
                     (lines[i].startswith(' ') or\
                     lines[i].startswith(field_identifier))):
                    # if continuation combine the lines while treating the
                    # spaces nicely, ie, multiple spaces -> one space
                    # this is mostly just important for the
                    # lines that are strings such as title
                    result = result.rstrip() + line_split + lines[i].lstrip()
                    i += 1
                break
            i += 1
        # return the field of interest   
        return result
        
class AAIndex1Parser(AAIndexParser):
    """ Parse AAIndex1 file & return it as dict of AAIndex1 objects"""

    def _parse_record(self, lines):
        """ Parse a single record and return it as a AAIndex1Record Object """
        # init all of the fields each time, this is so that
        # if fields are missing they don't get the value from the last
        # record
        id = None
        description = None
        LITDB = None
        authors = None
        title = None
        citations = None
        comments = None
        correlating = {}
        data = [None] * 20

        id = self._get_field('H', lines)
        description = self._get_field('D', lines)
        LITDB = self._get_field('R', lines)
        authors = self._get_field('A', lines)
        title = self._get_field('T', lines)
        citations = self._get_field('J', lines)
        comments = self._get_field('*', lines)
        correlating = self._parse_correlating(self._get_field('C', lines))
        data = self._parse_data(self._get_field('I', lines))

        return AAIndex1Record(id, description, LITDB, authors,\
                title, citations, comments, correlating, data)
                    

    def _parse_correlating(self, raw):
        """ Parse Correlating entries from the current record """
        keys = []
        values = []
        raw = raw.lstrip()
        # Split by white space
        data = re.split('\s*', raw)

        i=0
        while(i<len(data)):
            # If it's even it's a key
            if((i % 2) == 0):
                keys += [data[i]]
            # if it's not even it's a value
            else:
                # convert values to floats
                try:
                    values += [float(data[i])]
                except ValueError:
                    values += [data[i]]
            i += 1
        result = dict(zip(keys, values))
        return result

    def _parse_data(self, raw):
        """ Parse the data field from current record into a dict
        """
        # init for use in result
        keys = 'ARNDCQEGHILKMFPSTWYV'  
        values = []
        
        # get rid of leading white spaces, it makes../ the reg exp act weird
        raw = raw.lstrip()
        # split by any number/ types of white spaces
        data = re.split('\s*', raw)
        # convert the data to a float while checking for invlaid data,
        # specifically the string 'NA' is present sometimes instead of data
        for i in data:
            try:
                values += [float(i)]
            except ValueError:
                values += i

        result = dict(zip(keys, values))
        # return the dict
        return result


class AAIndex2Parser(AAIndexParser):
    """ Parse AAIndex2 file & return it as dict of AAIndex2 objects"""

    def _parse_record(self, lines):
        """ Parse a single record and return it as a AAIndex2Record Object """
        # Init attributes of each record each run through
       
        id = None
        description = None
        LITDB = None
        authors = None
        title = None
        citations = None
        comments = None
        rowscols = None
        data = []

        # Fill in the values
        id = self._get_field('H', lines)
        description = self._get_field('D', lines)
        LITDB = self._get_field('R', lines)
        authors = self._get_field('A', lines)
        title = self._get_field('T', lines)
        citations = self._get_field('J', lines)
        comments = self._get_field('*', lines)
        raw_data = self._get_field('M', lines)

        rowscols = self._parse_rowscols(raw_data[:raw_data.find('BRK')])
        try:
            data = self._parse_data(raw_data[raw_data.find('BRK')+3:],\
            rowscols[0], rowscols[1])
        except IndexError:
            return None

        return AAIndex2Record(id, description, LITDB, authors,\
            title, citations, comments, data)                       

    def _parse_data(self, raw, rows, cols):
        """ Parse the data field from current record into dict """
        # init result dict
        result = None
        # get rid of leading white spaces, it make the reg exp act weird
        raw = raw.lstrip()
        # split by any number/ types of white spaces
        data = re.split('\s*', raw)


        # If square matrix
        if len(data) == (len(rows)*len(cols)):
            result = dict.fromkeys(rows)
            i = 0
            for r in rows:
                new_row = dict.fromkeys(cols)
                for c in cols:
                    try:
                        new_row[c] = float(data[i])
                    except ValueError:
                        new_row[c] = data[i]
                    i+=1
                result[r] = new_row

        # else if LTM
        elif len(data) == (len(cols)+1) * len(rows)/2 :
            result = dict.fromkeys(rows)
            i = 0
            for r in rows:
                new_row = dict.fromkeys(cols)
                for c in cols:
                    if cols.find(c) <= rows.find(r):
                        try:
                            new_row[c] = float(data[i])
                        except ValueError:
                            new_row[c] = data[i]
                        i += 1
                result[r] = new_row                      
            
        return result

    def _parse_rowscols(self, raw):
        """ Returns two element list, 0: rows info, 1: cols info
        
            This parses the data out of the data description line
            for each record in AAIndex2 so we know what the data is that
            we are looking at.
        """
        p ='[rows|cols]\s=\s([^ \t\n\r\f\v,]*)'
        result = []
        result += re.findall(p, raw)
        return result


class AAIndexRecord(object):
    """ Abstract class, stores records from AAIndex files """

    def __init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data):
        """ Stores data for individual AAIndex entires """

        self.ID = str(id)
        self.Description = str(description)
        self.LITDBEntryNum = str(LITDB_entry_num)
        self.Authors = str(authors)
        self.Title = str(title)
        self.Citation = str(citation)
        self.Comments = str(comments)
        self.Data = data

    def _toSquareDistanceMatrix(self, include_stops=False):
        """ Converts AAIndex Data to square distance matrix

            This abstract method must be overwritten for each subclass.
            The interface must be identical across subclasses, must
            take self and return new square matrix (for now).
        """
        pass


    def toDistanceMatrix(self, include_stops=False):
        """ Builds a DistanceMatrix object based on self """
        data = self._toSquareDistanceMatrix(include_stops=include_stops)

        # If there is missing or invalid data, data will be None
        # if that's the case return None for easy detection, otherwise
        # return a new DistanceMatrix object
        if data:
            return DistanceMatrix(data=data, info=self)

        return None

class AAIndex1Record(AAIndexRecord):
    """ Stores records from AAIndex1, inherits from AAIndexRecord """

    def __init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments,
                  correlating, data):
        """ Stores data for individual AAIndex 1 entires """

        # Call init from super class
        AAIndexRecord.__init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data)

        self.Correlating = correlating

    def _toSquareDistanceMatrix(self, include_stops=False):
        """ AAIndex1 data to square distance matrix

        """
        keys = self.Data.keys()
        if include_stops : keys += '*'

        # build result dict top layer, start empty
        result = {}
        for r in keys:
            new_row = {}
            for c in keys:
                if (r == '*' or c == '*'):
                    new_row[c] = None
                else:
                    # Build the ditance matrix by subtracting the
                    # value of each aminoacid and then taking the
                    # absolute value.  If the data can not be
                    # turned into a float, it's not a number, so the data
                    # is invalid. Return None for easy detection
                    try:
                        new_row[c] =\
                            abs(float(self.Data[r])
                             - float(self.Data[c]))
                    except ValueError:
                        return None
            result[r] = new_row

        return result


class AAIndex2Record(AAIndexRecord):
    """ Stores records from AAIndex2, inherits from AAIndexRecord  """
    def __init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data):
        """ Stores data for individual AAIndex 2 entires """

        # Call init from super class
        AAIndexRecord.__init__(self, id,
                  description, LITDB_entry_num,
                  authors, title,
                  citation, comments, data)


    def _toSquareDistanceMatrix(self, include_stops=False):
        """ Returns data as a square matrix

            Note: This method is not currently functional,
            we are awaiting information on how to process data into
            a distance matrix

        """
        # create a new dict based on self.Data so we don't alter self.Data

        result = dict(self.Data)
        # Add in the new row of stop codon data
        if include_stops:
            stop_row = {}
            for i in result:
                stop_row.update({i:None})
            result.update({'*':stop_row})
            for i in result:
                result[i].update({'*':None})

        # Right now we are only dealing with square matrices
        return result

def AAIndexLookup(records):
    """ Build a dict of AAIndexObjects hashed by ID """    
    result = {}
    for r in records:
        result[r.ID] = r

    return result
        
def AAIndex1FromFiles(file):
    """ Taking a file or list of data return a dict of AAIndex1Objects """
    aap = AAIndex1Parser()
    return AAIndexLookup(aap(file))

def AAIndex2FromFiles(file):
    """ Taking a file or list of data return a dict of AAIndex2Objects """
    aap = AAIndex2Parser()
    return AAIndexLookup(aap(file))    

Woese_data = """//
H WOEC730101
D Polar requirement (Woese, 1973)
R PMID:4588588
A Woese, C.R.
T Evolution of genetic code
J Naturwiss. 60, 447-459 (1973)
C GRAR740102    0.960  HOPT810101    0.886  HOPA770101    0.876
  LEVM760101    0.872  PRAM900101    0.871  ROSM880101    0.844
  WOLS870101    0.841  KUHL950101    0.837  OOBM770103    0.835
  VINM940101    0.834  PARJ860101    0.821  FUKS010102    0.820
  FAUJ880110    0.812  OOBM770101    0.804  ROSM880102    0.801
  NADH010102   -0.800  CIDH920105   -0.800  MEIH800103   -0.802
  ISOY800102   -0.803  EISD860103   -0.803  ROSG850102   -0.804
  TANS770103   -0.806  RADA880101   -0.812  BIOV880102   -0.819
  WIMW960101   -0.821  NISK860101   -0.822  PONP800103   -0.823
  CIDH920104   -0.823  RADA880108   -0.825  BIOV880101   -0.829
  PONP800108   -0.831  SWER830101   -0.832  EISD860101   -0.838
  MAXF760102   -0.842  DESM900102   -0.847  FAUJ830101   -0.880
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
     7.0     9.1    10.0    13.0     5.5     8.6    12.5     7.9     8.4     4.9
     4.9    10.1     5.3     5.0     6.6     7.5     6.6     5.3     5.7     5.6
//
"""

def getWoeseDistanceMatrix():
    """ Return the Woese Polar Requirement Distance Matrix """
    aaindexObjects = AAIndex1FromFiles(Woese_data.split('\n'))
    distance_matrices = {}
    for m in aaindexObjects:
        distance_matrices[m] = aaindexObjects[m].toDistanceMatrix()

    return distance_matrices['WOEC730101']    


# _________________________________________________________________________________
_aaindex = dict()
_pymol_auto_arg_update = lambda: None


def search(pattern, searchtitle=True, casesensitive=False):
    '''
    Search for pattern in description and title (optional) of all records and
    return matched records as list. By default search case insensitive.
    '''
    whatcase = lambda i: i
    if not casesensitive:
        pattern = pattern.lower()
        whatcase = lambda i: i.lower()
    matches = []
    for record in _aaindex.itervalues():
        if pattern in whatcase(record.desc) or searchtitle and pattern in whatcase(record.title):
            matches.append(record)
    return matches


def grep(pattern):
    '''
    Search for pattern in title and description of all records (case
    insensitive) and print results on standard output.
    '''
    for record in search(pattern):
        print(record)


class Record:

    '''
    Amino acid index (AAindex) Record
    '''
    aakeys = 'ARNDCQEGHILKMFPSTWYV'

    def __init__(self):
        self.key = None
        self.desc = ''
        self.ref = ''
        self.authors = ''
        self.title = ''
        self.journal = ''
        self.correlated = dict()
        self.index = dict()
        self.comment = ''

    def extend(self, row):
        i = len(self.index)
        for x in row:
            self.index[self.aakeys[i]] = x
            i += 1

    def get(self, aai, aaj=None, d=None):
        assert aaj is None
        return self.index.get(aai, d)

    def __getitem__(self, aai):
        return self.get(aai)

    def median(self):
        x = sorted(filter(None, self.index.values()))
        half = len(x) / 2
        if len(x) % 2 == 1:
            return x[half]
        return (x[half - 1] + x[half]) / 2.0

    def __str__(self):
        desc = self.desc.replace('\n', ' ').strip()
        return '%s(%s: %s)' % (self.__class__.__name__, self.key, desc)


class MatrixRecord(Record):

    '''
    Matrix record for mutation matrices or pair-wise contact potentials
    '''

    def __init__(self):
        Record.__init__(self)
        self.index = []
        self.rows = dict()
        self.cols = dict()

    def extend(self, row):
        self.index.append(row)

    def _get(self, aai, aaj):
        i = self.rows[aai]
        j = self.cols[aaj]
        return self.index[i][j]

    def get(self, aai, aaj, d=None):
        try:
            return self._get(aai, aaj)
        except:
            pass
        try:
            return self._get(aaj, aai)
        except:
            return d

    def __getitem__(self, aaij):
        return self.get(aaij[0], aaij[1])

    def median(self):
        x = []
        for y in self.index:
            x.extend(filter(None, y))
        x.sort()
        if len(x) % 2 == 1:
            return x[len(x) / 2]
        return sum(x[len(x) / 2 - 1:len(x) / 2 + 1]) / 2.0


def get(key):
    '''
    Get record for key
    '''
    if len(_aaindex) == 0:
        init()
    return _aaindex[key]


def _float_or_None(x):
    if x == 'NA' or x == '-':
        return None
    return float(x)


def init(path=None, index='13'):
    '''
    Read in the aaindex files. You need to run this (once) before you can
    access any records. If the files are not within the current directory,
    you need to specify the correct directory path. By default all three
    aaindex files are read in.
    '''
    index = str(index)
    if path is None:
        for path in [os.path.split(__file__)[0], '.', cmd.get('fetch_path')]:
            if os.path.exists(os.path.join(path, 'aaindex' + index[0])):
                break
        print >> sys.stderr, 'path =', path
    if '1' in index:
        _parse(path + '/aaindex1', Record)
    if '2' in index:
        _parse(path + '/aaindex2', MatrixRecord)
    if '3' in index:
        _parse(path + '/aaindex3', MatrixRecord)
    _pymol_auto_arg_update()


def init_from_file(filename, type=Record):
    _parse(filename, type)


def _parse(filename, rec, quiet=True):
    '''
    Parse aaindex input file. `rec` must be `Record` for aaindex1 and
    `MarixRecord` for aaindex2 and aaindex3.
    '''
    if not os.path.exists(filename):
        import urllib
        url = 'ftp://ftp.genome.jp/pub/db/community/aaindex/' + os.path.split(filename)[1]
        print('Downloading "%s"' % (url))
        filename = urllib.urlretrieve(url, filename)[0]
        print('Saved to "%s"' % (filename))
    f = open(filename)

    current = rec()
    lastkey = None

    for line in f:
        key = line[0:2]
        if key[0] == ' ':
            key = lastkey

        if key == '//':
            _aaindex[current.key] = current
            current = rec()
        elif key == 'H ':
            current.key = line[2:].strip()
        elif key == 'R ':
            current.ref += line[2:]
        elif key == 'D ':
            current.desc += line[2:]
        elif key == 'A ':
            current.authors += line[2:]
        elif key == 'T ':
            current.title += line[2:]
        elif key == 'J ':
            current.journal += line[2:]
        elif key == '* ':
            current.comment += line[2:]
        elif key == 'C ':
            a = line[2:].split()
            for i in range(0, len(a), 2):
                current.correlated[a[i]] = float(a[i + 1])
        elif key == 'I ':
            a = line[1:].split()
            if a[0] != 'A/L':
                current.extend(map(_float_or_None, a))
            elif list(Record.aakeys) != [i[0] for i in a] + [i[-1] for i in a]:
                print('Warning: wrong amino acid sequence for', current.key)
            else:
                try:
                    assert list(Record.aakeys[:10]) == [i[0] for i in a]
                    assert list(Record.aakeys[10:]) == [i[2] for i in a]
                except:
                    print('Warning: wrong amino acid sequence for', current.key)
        elif key == 'M ':
            a = line[2:].split()
            if a[0] == 'rows':
                if a[4] == 'rows':
                    a.pop(4)
                assert a[3] == 'cols' and len(a) == 6
                i = 0
                for aa in a[2]:
                    current.rows[aa] = i
                    i += 1
                i = 0
                for aa in a[5]:
                    current.cols[aa] = i
                    i += 1
            else:
                current.extend(map(_float_or_None, a))
        elif not quiet:
            print('Warning: line starts with "%s"' % (key))

        lastkey = key

########## PYMOL ###########

# from Bio.SCOP.Raf import to_one_letter_code
# See also http://www.pymolwiki.org/index.php/Aa_codes
to_one_letter_code = {'PAQ': 'Y', 'AGM': 'R', 'ILE': 'I', 'PR3': 'C',
                      'GLN': 'Q', 'DVA': 'V', 'CCS': 'C', 'ACL': 'R', 'GLX': 'Z', 'GLY': 'G',
                      'GLZ': 'G', 'DTH': 'T', 'OAS': 'S', 'C6C': 'C', 'NEM': 'H', 'DLY': 'K',
                      'MIS': 'S', 'SMC': 'C', 'GLU': 'E', 'NEP': 'H', 'BCS': 'C', 'ASQ': 'D',
                      'ASP': 'D', 'SCY': 'C', 'SER': 'S', 'LYS': 'K', 'SAC': 'S', 'PRO': 'P',
                      'ASX': 'B', 'DGN': 'Q', 'DGL': 'E', 'MHS': 'H', 'ASB': 'D', 'ASA': 'D',
                      'NLE': 'L', 'DCY': 'C', 'ASK': 'D', 'GGL': 'E', 'STY': 'Y', 'SEL': 'S',
                      'CGU': 'E', 'ASN': 'N', 'ASL': 'D', 'LTR': 'W', 'DAR': 'R', 'VAL': 'V',
                      'CHG': 'A', 'TPO': 'T', 'CLE': 'L', 'GMA': 'E', 'HAC': 'A', 'AYA': 'A',
                      'THR': 'T', 'TIH': 'A', 'SVA': 'S', 'MVA': 'V', 'SAR': 'G', 'LYZ': 'K',
                      'BNN': 'A', '5HP': 'E', 'IIL': 'I', 'SHR': 'K', 'HAR': 'R', 'FME': 'M',
                      'PYX': 'C', 'ALO': 'T', 'PHI': 'F', 'ALM': 'A', 'PHL': 'F', 'MEN': 'N',
                      'TPQ': 'A', 'GSC': 'G', 'PHE': 'F', 'ALA': 'A', 'MAA': 'A', 'MET': 'M',
                      'UNK': 'X', 'LEU': 'L', 'ALY': 'K', 'SET': 'S', 'GL3': 'G', 'TRG': 'K',
                      'CXM': 'M', 'TYR': 'Y', 'SCS': 'C', 'DIL': 'I', 'TYQ': 'Y', '3AH': 'H',
                      'DPR': 'P', 'PRR': 'A', 'CME': 'C', 'IYR': 'Y', 'CY1': 'C', 'TYY': 'Y',
                      'HYP': 'P', 'DTY': 'Y', '2AS': 'D', 'DTR': 'W', 'FLA': 'A', 'DPN': 'F',
                      'DIV': 'V', 'PCA': 'E', 'MSE': 'M', 'MSA': 'G', 'AIB': 'A', 'CYS': 'C',
                      'NLP': 'L', 'CYQ': 'C', 'HIS': 'H', 'DLE': 'L', 'CEA': 'C', 'DAL': 'A',
                      'LLP': 'K', 'DAH': 'F', 'HMR': 'R', 'TRO': 'W', 'HIC': 'H', 'CYG': 'C',
                      'BMT': 'T', 'DAS': 'D', 'TYB': 'Y', 'BUC': 'C', 'PEC': 'C', 'BUG': 'L',
                      'CYM': 'C', 'NLN': 'L', 'CY3': 'C', 'HIP': 'H', 'CSO': 'C', 'TPL': 'W',
                      'LYM': 'K', 'DHI': 'H', 'MLE': 'L', 'CSD': 'A', 'HPQ': 'F', 'MPQ': 'G',
                      'LLY': 'K', 'DHA': 'A', 'DSN': 'S', 'SOC': 'C', 'CSX': 'C', 'OMT': 'M',
                      'DSP': 'D', 'PTR': 'Y', 'TRP': 'W', 'CSW': 'C', 'EFC': 'C', 'CSP': 'C',
                      'CSS': 'C', 'SCH': 'C', 'OCS': 'C', 'NMC': 'G', 'SEP': 'S', 'BHD': 'D',
                      'KCX': 'K', 'SHC': 'C', 'C5C': 'C', 'HTR': 'W', 'ARG': 'R', 'TYS': 'Y',
                      'ARM': 'R', 'DNP': 'A'}


def aaindex2b(key='KYTJ820101', selection='(all)', quiet=0, var='b'):
    '''
DESCRIPTION

    "aaindex" looks up the Amino Acid Index from
      http://www.genome.jp/aaindex/
    for the given key and assignes b-factors to the given selection. Unknown
    residues get the average index value assigned.

USAGE

    aaindex2b [key [, selection]]

ARGUMENTS

    key = string: Key of AAindex entry

    selection = string: atoms to assign b-factors {default: (all)}

EXAMPLE

    # Hydropathy index by Kyte-Doolittle
    aaindex2b KYTJ820101
    spectrumany b, white yellow forest
    show surface
    '''
    from pymol import cmd, stored
    entry = get(key)
    median = entry.median()

    if int(quiet) != 0:
        print(entry.desc.strip())

    def lookup(resn):
        one_letter = to_one_letter_code.get(resn, 'X')
        value = entry.get(one_letter)
        if value is None:
            return median
        return value
    stored.aaindex = lookup

    cmd.alter(selection, var + '=stored.aaindex(resn)')


def pmf(key, cutoff=7.0, selection1='(name CB)', selection2='', state=1, quiet=1):
    '''
DESCRIPTION

    Potential of Mean Force

ARGUMENTS

    key = string: aaindex key

    cutoff = float: distance cutoff {default: 7.0}
    cutoff = (float, float): distance shell

    selection1 = string: atom selection {default: (name CB)}

    selection2 = string: atom selection {default: selection1}

NOTES

    Does also support a list of keys and a list of cutoffs to deal with
    multiple distance shells.

EXAMPLES

    # distance dependent c-beta contact potentials
    pmf SIMK990101, 5,         /2x19//A//CB
    pmf SIMK990102, [5, 7.5],  /2x19//A//CB
    pmf [SIMK990101, SIMK990102, SIMK990103], [0, 5, 7.5, 10], /2x19//A//CB

    # interface potential
    sidechaincenters 2x19_scc, 2x19
    pmf KESO980102, 7.0, /2x19_scc//A, /2x19_scc//B
    distance /2x19_scc//A, /2x19_scc//B, cutoff=7.0
    '''
    from pymol import cmd, stored
    from chempy import cpv
    if cmd.is_string(key):
        if key.lstrip().startswith('['):
            key = cmd.safe_alpha_list_eval(key)
        else:
            key = [key]
    if cmd.is_string(cutoff):
        cutoff = eval(cutoff)
    if not cmd.is_sequence(cutoff):
        cutoff = [cutoff]
    if len(cutoff) == len(key):
        cutoff = [0.0] + list(cutoff)
    if len(cutoff) != len(key) + 1:
        print('Error: Number of keys and number of cutoffs inconsistent')
        return
    state = int(state)
    quiet = int(quiet)
    if len(selection2) == 0:
        selection2 = selection1
    if not quiet and len(key) > 1:
        print('Distance shells:')
        for i in range(len(key)):
            print('%s %.1f-%.1f' % (key[i], cutoff[i], cutoff[i + 1]))

    idmap = dict()
    cmd.iterate_state(state, '(%s) or (%s)' % (selection1, selection2),
                      'idmap[model,index] = [(resn,name),(x,y,z)]', space={'idmap': idmap})
    twoN = cmd.count_atoms(selection1) + cmd.count_atoms(selection2)
    pairs = cmd.find_pairs(selection1, selection2, cutoff=max(cutoff),
                           state1=state, state2=state)
    if len(pairs) == 0:
        print('Empty pair list')
        return 0.0

    matrix = map(get, key)
    for i in matrix:
        assert isinstance(i, MatrixRecord)

    i_list = range(len(key))
    u_sum = 0
    count = 0
    for id1, id2 in pairs:
        a1 = idmap[id1]
        a2 = idmap[id2]
        r = cpv.distance(a1[1], a2[1])
        for i in i_list:
            if cutoff[i] <= r and r < cutoff[i + 1]:
                try:
                    aa1 = to_one_letter_code[a1[0][0]]
                    aa2 = to_one_letter_code[a2[0][0]]
                    u_sum += matrix[i].get(aa1, aa2)
                    count += 1
                except:
                    print('Failed for', a1[0], a2[0])

    value = float(u_sum) / twoN
    if not quiet:
        print('PMF: %.4f (%d contacts, %d residues)' % (value, count, twoN))
    return value

try:
    from pymol import cmd
    cmd.extend('aaindex2b', aaindex2b)
    cmd.extend('pmf', pmf)

    def pymol_auto_arg_update():
        aaindexkey_sc = cmd.Shortcut(_aaindex.keys())
        cmd.auto_arg[0].update({
            'aaindex2b': [aaindexkey_sc, 'aaindexkey', ', '],
            'pmf': [aaindexkey_sc, 'aaindexkey', ', '],
        })
        cmd.auto_arg[1].update({
            'aaindex2b': [cmd.selection_sc, 'selection', ''],
        })
        cmd.auto_arg[2].update({
            'pmf': [cmd.selection_sc, 'selection', ''],
        })
        cmd.auto_arg[3].update({
            'pmf': [cmd.selection_sc, 'selection', ''],
        })
    _pymol_auto_arg_update = pymol_auto_arg_update
except:
    pass

# vi: ts=4:sw=4:smarttab:expandtab

def importRecord(entry_id):
    # return an entry using accession number
    
    """
    usage example:
       aaIndex, key, desc, ref, authors, title, journal, correlated, index, comment = importRecord('KRIW790103')
       print aaIndex
       print key
       print desc
       print ref
       print authors
    """
    #aaindex.init(path='../aaindex')  #database files
    init(path='../aaindex')  #database files
    aakeys = 'ARNDCQEGHILKMFPSTWYV'
    #aaindex.grep("")
    #x = aaindex.get(entry_id)
    x = get(entry_id)
    aaIndex={}
    for aa in aakeys:
       aaIndex[aa] = x.get(aa)
    #key = x.key
    #desc = x.desc
    #ref = x.ref
    #authors = x.authors
    #title = x.title
    #journal = x.journal
    #correlated = x.correlated #dict()
    #index = x.index           #dict()
    #comment = x.comment 
    return aaIndex #, key, desc, ref, authors, title, journal, correlated, index, comment   

def getRecord(entry_id):
    #aaindex.init(path='../aaindex')
    init(path='../aaindex')
    #from this you can get: key, desc, ref, authors, title, journal, correlated, index, comment
    #e.g. record.key , record.desc ...etc
    record = get(entry_id)
    return record


def importEntryIDs():
    #import all accession numbers from the aaindex file as a list
    """
    usage example:
       nums = importEntryIDs()
       for i in nums:
          print i
       print len(nums)
    """ 
    f = open('../aaindex/aaindex1')
    accessionNumbers = list(AAIndex1FromFiles(f))
    return accessionNumbers

    

def test():
       record = getRecord('KRIW790103')
       print("key: " + record.key)
       print("desc: " + record.desc)
       print("ref: " + record.ref)
       print("author: " + record.authors)
#if __name__ == "__main__":
   #test()


###############################################################################
AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

#c Cys is separated from class 3 because of its ability to form disulfide bonds. 

_repmat={1:["A",'G','V'],2:['I','L','F','P'],3:['Y','M','T','S'],4:['H','N','Q','W'],5:['R','K'],6:['D','E'],7:['C']}

###############################################################################


class Proteindescriptors(object):
   def __init__(self):
       super(Proteindescriptors, self).__init__()   


   #****** Getting properties  *****

   def getAAIndex(self, ID):
      '''
       Given AAindex ID will generate the index as a dictionary
      '''
      from propy import AAIndex
      return AAIndex.GetAAIndex1(ID)


   def getHydrophobicity(self):
       '''
        To get Hydrophobicity
       '''
       from propy import PseudoAAC
       Hydrophobicity = PseudoAAC._Hydrophobicity
       return Hydrophobicity

   def getpK1(self):
       '''
        To get pK1
       '''
       from propy import PseudoAAC
       pK1 = PseudoAAC._pK1
       return pK1

   def getpK2(self):
       '''
        To get pK2
       '''
       from propy import PseudoAAC
       pK2 = PseudoAAC._pK2
       return pK2


   def getResidueMass(self):
       '''
        To get mass
       '''
       from propy import PseudoAAC
       residueMass = PseudoAAC._residuemass
       return residueMass

   def getHydrophilicity(self):
       '''
        Hydrophilicity
       '''
       from propy import PseudoAAC
       hydrophilicity = PseudoAAC._hydrophilicity
       return hydrophilicity

   def getPI(self):
       '''
        Isoelectric point
       '''
       from propy import PseudoAAC
       pI = PseudoAAC._pI
       return pI

   #******************* Conjoint triad *****************************
   def _Str2Num(self, proteinsequence):
        """
        translate the amino acid letter into the corresponding class based on the
        
        given form.
        
        """
        repmat={}
        for i in _repmat:
            for j in _repmat[i]:
                repmat[j]=i

        res=proteinsequence
        for i in repmat:
            res=res.replace(i,str(repmat[i]))
        return res


   ###############################################################################
   def CalculateConjointTriad(self, proteinsequence):

        """
        Calculate the conjoint triad features from protein sequence.
        
        Useage:
        
        res = CalculateConjointTriad(protein)
        
        Input: protein is a pure protein sequence.
        
        Output is a dict form containing all 343 conjoint triad features.
        """
        res={}
        proteinnum=self._Str2Num(proteinsequence)
        for i in range(8):
                for j in range(8):
                        for k in range(8):
                                temp=str(i)+str(j)+str(k)
                                res[temp]=proteinnum.count(temp)
        return res


   #**********************************Shannon entropy ***************
   def EntropyJunk(self, seq):
        import math
        log2=lambda x:math.log(x)/math.log(2)
        exr={}
        infoc=0.0
        for each in seq:
            try:
                exr[each]+=1
            except:
                exr[each]=1
        seqlen=len(seq)
        for k,v in exr.items():
            freq  =  1.0*v/seqlen
            infoc+=freq*log2(freq)
        infoc*=-1
        return infoc


   def calRelativeEntropyJunk(self, seq):
         from math import log
         resCodes = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
         #resCodes = "GCAT"
         N = float(len(seq))
         base = 1.0 / len(resCodes)
         prop ={}
         for r in resCodes:
            prop[r] = 0
         for r in seq:
            prop[r]+=1
         for r in resCodes:
            prop[r] /= N

         H = 0
         for r in resCodes:
             if prop[r] != 0.0:
                 h = prop[r] * log(prop[r] / base, 2.0)
                 H += h
         H /= log(base, 2.0)
         return H   

   def CalInformationGainJunk(self, seq):
       IG = self.Entropy(seq) - self.calRelativeEntropy(seq)

   def getCKSAAP(self, proseq):
       AALetter2 = AALetter
       res = {}
       AANum = 0
       for k in range(0,6):
           K = "X"*k
           for i, AA in enumerate(AALetter):
               for j, AA2 in enumerate(AALetter2):
                   searchPattern = AA+'[a-zA-Z]{'+str(k)+'}'+AA2
                   AANum = AANum + 1
                   resValue = len(re.findall(searchPattern,proseq))
                   resKey = AA+K+AA2+"_"+str(AANum)
                   res[resKey] = resValue
       return(res)

   #************** Discrete model: Ammino Acid Composition **********
   def getAAComposition(self, proseq):
       '''
       * Amino acid composition
       * 20 feature / sequence
       * The relative frequency is multiplied by 100
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetAAComp()

   def getDipeptideComp(self, proseq):
       '''
       * Dipeptide composition
       * 400 feature / sequence
       '''
       from propy import AAComposition as AAC
       return AAC.CalculateDipeptideComposition(proseq)
       
   def getAADipeptideComposition(self, proseq):
       '''
       * AADs, dipeptide and 3-mers
       *  8420 Features / sequence
       '''
       from propy import AAComposition as AAC
       return AAC.CalculateAADipeptideComposition(proseq)

   def getSpectrumDict(self, proseq):
       '''
       * Calcualte the spectrum descriptors of 3-mers for a given protein
       *  8000 Features / sequence
       '''
       from propy import AAComposition as AAC
       return AAC.GetSpectrumDict(proseq)

   #************** End of Discrete model: Ammino Acid Composition **********

   #****** Automated active site detection, docking, and scoring (AADs)  *****

   def getCTD(self, proseq):
       '''
       * CTD: Composition Translation Distribution
       * Calculate all CTD descriptors based seven different properties of AADs. 
       * 21 feature / sequence
       * Composition
       '''
       from propy import CTD
       return CTD.CalculateCTD(proseq)

   def getComposition(self, proseq):
       '''
       * Calculate all composition descriptors based seven different properties of
       *  automated active site detection, docking, and scoring (AADs) 
       * 21 feature / sequence
       * Composition
       '''
       from propy import CTD
       return CTD.CalculateC(proseq)
   
   def getTranslation(self, proseq):
       '''
       * Calculate all transition descriptors based seven different properties of AADs 
       * 21 feature / sequence
       * Translation
       '''
       from propy import CTD
       return CTD.CalculateT(proseq)

   def getDistribution(self, proseq):
       '''
       * Calculate all distribution descriptors based seven different properties of AADs 
       * 21 feature / sequence
       * Distribution
       '''
       from propy import CTD
       return CTD.CalculateD(proseq)

   def getCalculateCompositionCharge(self, proseq):
       '''
       * calculating composition descriptors based on Charge of AADs 
       * 3 feature / sequence      
       '''
       from propy import CTD
       return CTD.CalculateCompositionCharge(proseq)


   def getCalculateCompositionHydrophobicity(self, proseq):
       '''
       * calculating composition descriptors based on Hydrophobicity of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateCompositionHydrophobicity(proseq)

   def getCalculateCompositionNormalizedVDWV(self, proseq):
       '''
       * composition descriptors based on NormalizedVDWV of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateCompositionNormalizedVDWV(proseq)


   def getCalculateCompositionPolarity(self, proseq):
       '''
       * calculating composition descriptors based on Polarity of of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateCompositionPolarity(proseq)


   def getCalculateCompositionPolarizability(self, proseq):
       '''
       * calculating composition descriptors based on Polarizability of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateCompositionPolarizability(proseq)


   def getCalculateCompositionSecondaryStr(self, proseq):
       '''
       * calculating composition descriptors based on SecondaryStr of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateCompositionSecondaryStr(proseq)

   def getCalculateCompositionSolventAccessibility(self, proseq):
       '''
       * calculating composition descriptors based on SolventAccessibilit of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateCompositionSolventAccessibility(proseq)


   def getCalculateDistributionCharge(self, proseq):
       '''
       * calculating Distribution descriptors based on Charge of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionCharge(proseq)


   def getCalculateDistributionHydrophobicity(self, proseq):
       '''
       * calculating Distribution descriptors based on Hydrophobicity of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionHydrophobicity(proseq)


   def getCalculateDistributionNormalizedVDWV(self, proseq):
       '''
       * calculating Distribution descriptors based on NormalizedVDWV of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionNormalizedVDWV(proseq)

   def getCalculateDistributionPolarity(self, proseq):
       '''
       * calculating Distribution descriptors based on Polarity of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionPolarity(proseq)


   def getCalculateDistributionPolarizability(self, proseq):
       '''
       * calculating Distribution descriptors based on Polarizability of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionPolarizability(proseq)


   def getCalculateDistributionSecondaryStr(self, proseq):
       '''
       * calculating Distribution descriptors based on SecondaryStr of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionSecondaryStr(proseq)


   def getCalculateDistributionSolventAccessibility(self, proseq):
       '''
       * calculating Distribution descriptors based on SolventAccessibility of AADs 
       * 15 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateDistributionSolventAccessibility(proseq)


   def getCalculateTransitionCharge(self, proseq):
       '''
       * calculating Transition descriptors based on Charge of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionCharge(proseq)


   def getCalculateTransitionHydrophobicity(self, proseq):
       '''
       * calculating Transition descriptors based on Hydrophobicity of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionHydrophobicity(proseq)


   def getCalculateTransitionNormalizedVDWV(self, proseq):
       '''
       * calculating Transition descriptors based on NormalizedVDWV of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionNormalizedVDWV(proseq)


   def getCalculateTransitionPolarity(self, proseq):
       '''
       * calculating Transition descriptors based on Polarity of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionPolarity(proseq)


   def getCalculateTransitionPolarizability(self, proseq):
       '''
       * calculating Transition descriptors based on Polarizability of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionPolarizability(proseq)


   def getCalculateTransitionSecondaryStr(self, proseq):
       '''
       * calculating Transition descriptors based on SecondaryStr of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionSecondaryStr(proseq)


   def getCalculateTransitionSolventAccessibility(self, proseq):
       '''
       * calculating Transition descriptors based on SolventAccessibility of AADs 
       * 3 feature / sequence     
       '''
       from propy import CTD
       return CTD.CalculateTransitionSolventAccessibility(proseq)


   #****** End of Automated active site detection, docking, and scoring (AADs)  *****


   #****** Correlation and pseudo amino acids  *****

   def getGearyAutoCorr(self, proseq):
       '''
        Geary autocoorrelation
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetGearyAuto()

   def getcusGearyAutop1(self, proseq, ID):
       '''
        Custom Geary autocoorrelation
       '''
       from propy.PyPro import GetProDes
       import propy
       from propy import AAIndex
       #proindex = AAIndex.GetAAIndex1(ID, path=propy.__path__[0])
       proindex = AAIndex.GetAAIndex1(ID)
       Des = GetProDes(proseq)
       return Des.GetGearyAutop(proindex)
     
   def getcusGearyAutop2(self, proseq, AAP):
       '''
        Custom Geary autocoorrelation
        APP is a dict
       '''
       from propy.PyPro import GetProDes
       import propy
       from propy import AAIndex
       Des = GetProDes(proseq)
       return Des.GetGearyAutop(AAP=AAP)


   def getGetMoranAuto(self, proseq):
       '''
        Moran autocorrelation descriptors (240)
        Features: 240 / protein
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetMoranAuto()

   def getGetMoranAutop1(self, proseq, ID):  
       '''
        Moran autocorrelation descriptors for the given property
        Features: 30 / protein
       '''
       from propy.PyPro import GetProDes
       import propy
       from propy import AAIndex
       proindex = AAIndex.GetAAIndex1(ID, path=propy.__path__[0])
       Des = GetProDes(proseq)
       return Des.GetMoranAutop(proindex)


   def getGetMoranAutop2(self, proseq, AAP):  
       '''
        Moran autocorrelation descriptors for the given property
        Features: 30 / protein
        AAP is a dict e.g. Hydrophobicity
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetMoranAutop(AAP=AAP)


   def getGetMoreauBrotoAuto(self, proseq):
       '''
        Normalized Moreau-Broto autocorrelation descriptors (240)
        Features: 240 / protein
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetMoreauBrotoAuto()

   def getGetMoreauBrotoAutop1(self, proseq, ID):  
       '''
        Normalized Moreau-Broto autocorrelation descriptors for the given property (30)
        Features: 30 / protein
       '''
       from propy.PyPro import GetProDes
       import propy
       from propy import AAIndex
       proindex = AAIndex.GetAAIndex1(ID, path=propy.__path__[0])
       Des = GetProDes(proseq)
       return Des.GetMoreauBrotoAutop(proindex)


   def getGetMoreauBrotoAutop2(self, proseq, AAP):  
       '''
        Normalized Moreau-Broto autocorrelation descriptors for the given property (30)
        Features: 30 / protein
        AAP is a dict e.g. Hydrophobicity
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetMoreauBrotoAutop(AAP=AAP)



   def getGetPAAC(self, proseq, lamda, weight):
       '''
        Type I Pseudo amino acid composition descriptors (default is 30)
        lamda < protein ength
        weight => 0.05 - 0.7
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetPAAC(lamda=lamda, weight=weight)


   def getGetPAACp1(self, proseq, ID, lamda, weight):
       '''
        Custom Pseudo amino acid composition
        ID is AAIndex ID
       '''
       from propy.PyPro import GetProDes
       import propy
       from propy import AAIndex
       proindex = AAIndex.GetAAIndex1(ID, path=propy.__path__[0])
       Des = GetProDes(proseq)
       return Des.GetPAACp(AAP=[proindex], lamda=lamda, weight=weight)

   def getGetPAACp2(self, proseq, AAP, lamda, weight):
       '''
        Type I Pseudo amino acid composition descriptors for the given properties (default is 30)
        AAP is a list
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetPAACp(AAP=AAP, lamda=lamda, weight=weight)

   def getGetAPseudoAAC1(self, proseq, lamda, weight):
       '''
        Computing the first 20 of type II pseudo-amino acid compostion descriptors based on
        [_Hydrophobicity,_hydrophilicity].

       '''
       from propy import PseudoAAC
       return PseudoAAC.GetAPseudoAAC1(proseq, lamda=lamda, weight=weight)
       
   #**********************************************************************
   def getGetSequenceOrderCorrelationFactorOLD(self, proseq, k=1, AAP=[]):
       '''
        Computing the Sequence order correlation factor with gap equal to k based on 
        the given properities.
        k is the gap
       '''
       from propy import PseudoAAC as AAC
       return AAC.GetSequenceOrderCorrelationFactor(proseq, k=k, AAP=AAP)

   def getGetSequenceOrderCorrelationFactor(self, proseq, ID1, ID2, k=1):
       '''
        Computing the Sequence order correlation factor with gap equal to k based on 
        the given properities.
        k is the gap
       '''
       from propy import PseudoAAC as AAC
       import propy
       from propy import AAIndex
       proindex1 = AAIndex.GetAAIndex1(ID1, path=propy.__path__[0])
       proindex2 = AAIndex.GetAAIndex1(ID2, path=propy.__path__[0])
       #proindex1 = self.getAAIndex(ID1)
       #proindex2 = self.getAAIndex(ID2)
       AAP = []
       AAP.append(proindex1)
       AAP.append(proindex2)
       return AAC.GetSequenceOrderCorrelationFactor(proseq, k=k, AAP=AAP)
   #*********************************************************************
   def distMatrix(self, csvFile):
       aa = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
       matrixFile = csvFile
       f = open(matrixFile, "r")
       i = 0
       j = 0
       mat = {}
       for line in f.readlines():
          line = line.split(",")
          for j in range(20):
             key = aa[i] + aa[j]
             val = float(line[j])
             j = j +1
             mat[key] = val
          i = i +1
       return mat

       
   def getGetSequenceOrderCorrelationFactorForAPAAC(self, proseq, k):
       '''
        Computing the Sequence order correlation factor with gap equal to k based on 
        [_Hydrophobicity,_hydrophilicity] for APAAC (type II PseAAC)
        k is the gap
       '''
       from propy import PseudoAAC as AAC
       return AAC.GetSequenceOrderCorrelationFactorForAPAAC(proseq, k=k)




   def getGetQSO(self, proseq, maxlag, weight):
       '''
        Quasi sequence order descriptors  default is 50
        maxlag is the maximum lag and the length of the protein should be larger than maxlag. default is 45
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetQSO(maxlag=maxlag, weight=weight)

   def getGetQSOp(self, proseq, distancematrix, maxlag, weight):
       '''
        Quasi sequence order descriptors  default is 50
        maxlag is the maximum lag and the length of the protein should be larger than maxlag. default is 45
        distancematrix is a dict form containing 400 distance values
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetQSOp(distancematrix=distancematrix, maxlag=maxlag, weight=weight)

   def getGetSOCN(self, proseq, maxlag):
       '''
        Sequence order coupling numbers  default is 45
        maxlag is the maximum lag and the length of the protein should be larger than maxlag. default is 45.
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetSOCN(maxlag=maxlag)

   def getGetSOCNp(self, proseq, distancematrix, maxlag):
       '''
        Quasi sequence order descriptors  default is 50
        maxlag is the maximum lag and the length of the protein should be larger than maxlag. default is 45
        distancematrix is a dict form containing 400 distance values
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetSOCNp(distancematrix=distancematrix, maxlag=maxlag)

   # --------------------- those return string -----------------------------
   def getGetSubSeq(self, proseq, ToAA='S', window=3):
       '''
        obtain the sub sequences wit length 2*window+1, whose central point is ToAA
        ToAA is the central (query point) amino acid in the sub-sequence.
        window is the span.
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetSubSeq(ToAA=ToAA, window=window)


   def getGetTPComp(self, proseq):
       '''
        tri-peptide composition descriptors (8000)
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetTPComp()
   # --------------------- -----------------------------------------------------

   def getALLDesc(self, proseq):
       '''
        Get all protein descriptors
       '''
       from propy.PyPro import GetProDes
       Des = GetProDes(proseq)
       return Des.GetALL()


   #****** End Correlation and pseudo amino acids  *****


class PDescriptorsProcessing(object):
   def __init__(self):
       super(PDescriptorsProcessing, self).__init__()   

   def countFastaSeq(self, fastaFile):
        """
        Purpose: Count number of sequence in fasta multiple sequence file
        Argument: fasta file 
        Return: integer as number of sequences
        """
        count = 0
        '''
        #for line in file(fastaFile, 'r'):
        '''
        for line in open(fastaFile, 'r'):
           if line.startswith('>'):
              count +=1
        return count 


   def convertFastasToCsv0(self, fun, fastaFilespath, tempPath, fastaFiles):
         """The function has only one argument which is the seq """
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os 
         PAAC = {}                
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq) 
                      if temp > 0:
                         PAAC = fun(seq)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                      f2.write(PACCstr)


             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            line = line.rstrip('\n')
                            line += f2.readline()
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(PAA + '\n')
             f3.close()    

   def convertFastasToCsvID(self, fun, fastaFilespath, tempPath, fastaFiles, ID):
         """The function has only two arguments: seq and ID"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             #f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, ID)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                         #PACCstr = str(PAAC.values()).strip('[]') + ','
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, ID)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                      #PACCstr = str(PAAC.values()).strip('[]') + ','
                      f2.write(PACCstr)
             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            line = line.rstrip('\n')
                            line += f2.readline()
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(ID+'_'+PAA + '\n')
             f3.close()    



   def convertFastasToCsv2IDs(self, fun, fastaFilespath, tempPath, fastaFiles, ID1, ID2):
         """The function has only two arguments: seq and AAP"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         PAAC ={}
         nFile = len(fastaFiles)
         for i in range(nFile):
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)
                      if temp > 0:
                         PAAC = fun(seq, ID1, ID2)
                         PACCstr = str(PAAC) + ',\n'
                         f2.write(PACCstr)

             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)
                  if temp > 0:
                      PAAC = fun(seq, ID1, ID2)
                      PACCstr = str(PAAC) + ',\n'
                      f2.write(PACCstr)


             f2.close()
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            line = line.rstrip('\n')
                            line += f2.readline()
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(ID1+'_'+ID2+'_'+PAA + '\n')
             f3.close()


   def convertFastasToCsvAAP(self, fun, fastaFilespath, tempPath, fastaFiles, AAP, pname):
         """The function has only two arguments: seq and AAP"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         #file_to_delete = tempPath + "*.fastaHA"
         #if os.path.exists(file_to_delete):
             #os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] +"HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, AAP)
                         #PACCstr = str(PAAC.values()).strip('[]') + ','
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, AAP)
                      #PACCstr = str(PAAC.values()).strip('[]') + ','
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                      f2.write(PACCstr)


             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            line = line.rstrip('\n')
                            line += f2.readline()
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(pname+'_'+PAA + '\n')
             f3.close()    

   def convertFastasToCsvLamdaWeight(self, fun, fastaFilespath, tempPath, fastaFiles, lamda, weight):
         """The function has only two arguments: seq and lamda, weight"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         des =[]
         PAAC = {}
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, lamda=lamda, weight=weight)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                         f2.write(PACCstr)
                         des = PAAC.keys()
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, lamda=lamda, weight=weight)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                      f2.write(PACCstr)
                      des = PAAC.keys()
                      
             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                           
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             #des = PAAC.keys()
             for PAA in des:
                   f3.write(PAA + '\n')
             f3.close()    
   def convertFastasToCsvLamdaWeightID(self, fun, fastaFilespath, tempPath, fastaFiles, lamda, weight, ID):
         """The function has only two arguments: seq and lamda, weight"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         des = []
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, ID, lamda=lamda, weight=weight)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                         f2.write(PACCstr)
                         des = PAAC.keys()
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, ID, lamda=lamda, weight=weight)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                      f2.write(PACCstr)
                      des = PAAC.keys()

             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             #des = PAAC.keys()
             for PAA in des:
                   f3.write(ID+'_'+PAA + '\n')
             f3.close()    

   def convertFastasToCsvLamdaWeightAAP(self, fun, fastaFilespath, tempPath, fastaFiles, lamda, weight, AAP, pname):
         """The function has only two arguments: seq and lamda, weight"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, AAP=AAP, lamda=lamda, weight=weight)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, AAP, lamda=lamda, weight=weight)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                      f2.write(PACCstr)


             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(pname+'_'+PAA + '\n')
             f3.close()    



   def convertFastasToCsvWeightMaxlag(self, fun, fastaFilespath, tempPath, fastaFiles, weight, maxlag):
         """The function has only two arguments: seq and weight, length"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, weight=weight, maxlag=maxlag)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, weight=weight, maxlag=maxlag)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                      f2.write(PACCstr)


             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(PAA + '\n')
             f3.close()    

   def convertFastasToCsvWeightMaxlagMatrix(self, fun, fastaFilespath, tempPath, fastaFiles, weight, maxlag, matrix, mat):
         """The function has only two arguments: seq and weight, length"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         PAAC = {}
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         nFile = len(fastaFiles)
         for i in range(nFile):
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)
                      if temp > 0:
                         PAAC = fun(seq, matrix, maxlag=maxlag, weight=weight)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)
                  if temp > 0:
                      PAAC = fun(seq, matrix, maxlag=maxlag, weight=weight)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                      f2.write(PACCstr)


             f2.close()
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(mat+PAA + '\n')
             f3.close()




   def convertFastasToCsvMaxlag(self, fun, fastaFilespath, tempPath, fastaFiles, maxlag):
         """The function has only two arguments: seq and length"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, maxlag=maxlag)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, maxlag=maxlag)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                      f2.write(PACCstr)


             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(PAA + '\n')
             f3.close()    


   def convertFastasToCsvMaxlagMatrix(self, fun, fastaFilespath, tempPath, fastaFiles, maxlag, matrix, mat):
         """The function has only two arguments: seq and length"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         PAAC = {}
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)

         nFile = len(fastaFiles)
         for i in range(nFile):
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)
                      if temp > 0:
                         PAAC = fun(seq, matrix, maxlag=maxlag)
                         PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                         f2.write(PACCstr)
             elif seqCount == 1:

                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)
                  if temp > 0:
                      PAAC = fun(seq, matrix, maxlag=maxlag)
                      PACCstr = str(PAAC.values()).strip('[]') + ',\n'
                      f2.write(PACCstr)

             f2.close()
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = PAAC.keys()
             for PAA in des:
                   f3.write(mat+PAA + '\n')
             f3.close()



   def convertFastasToCsvkspace(self, fun, fastaFilespath, tempPath, fastaFiles, kspace):
         """The function has only two arguments: seq and space"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PAAC = fun(seq, kspace)
                         PACCstr = str(PAAC).strip('[]') + ',\n' 
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                      PAAC = fun(seq, kspace)
                      PACCstr = str(PAAC).strip('[]') + ',\n' 
                      f2.write(PACCstr)


             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             #des = PAAC.keys()
             for i in range(len(PAAC)):
                   f3.write('CF' + str(i+1) + '\n')
             f3.close()    

   def convertFastasToCsvkspaceAAP(self, fun, fastaFilespath, tempPath, fastaFiles, lamda, AAP):
         """The function has only two arguments: seq and space
             fun output is not a dict it is a list"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         PACCstr=""
                         for kspace in range(1, lamda+1):
                            PAAC = fun(seq, kspace, AAP)
                            PACCstr = PACCstr + str(PAAC).strip('[]') + ','
                         PACCstr = PACCstr + '\n' 
                         f2.write(PACCstr)
             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                         PACCstr=""
                         for kspace in range(1, lamda+1):
                            PAAC = fun(seq, kspace, AAP)
                            PACCstr = PACCstr + str(PAAC).strip('[]') + ','
                         PACCstr = PACCstr + '\n' 
                         f2.write(PACCstr)
             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             for i in range(1, lamda+1):
                   f3.write('CF' + str(i) + '\n')
             f3.close()    

   #****************************************** Entropy zone ************************
   def Entropy(self, seq):
        import math
        log2=lambda x:math.log(x)/math.log(2)
        exr={}
        infoc=0
        for each in seq:
            try:
                exr[each]+=1
            except:
                exr[each]=1
        seqlen=len(seq)
        for k,v in exr.items():
            freq  =  1.0*v/seqlen
            infoc+=freq*log2(freq)
        infoc*=-1
        return infoc


   def calRelativeEntropy(self, seq):
         from math import log
         resCodes = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
         #resCodes = "GCAT"
         N = float(len(seq))
         base = 1.0 / len(resCodes)
         prop ={}
         for r in resCodes:
            prop[r] = 0
         for r in seq:
            prop[r]+=1
         for r in resCodes:
            prop[r] /= N

         H = 0
         for r in resCodes:
             if prop[r] != 0.0:
                 h = prop[r] * log(prop[r] / base, 2.0)
                 H += h
         #H /= log(base, 2.0)
         return H   


   def CalInformationGain(self, seq):
         IG = self.Entropy(seq) - self.calRelativeEntropy(seq)



   def convertFastasToCsvEntropy_relEntropy(self, fastaFilespath, tempPath, fastaFiles, FLAG):
         """The function has only two arguments: seq and weight, length"""
         #fun is a function that return a descriptor dict
         #This function convert each fasta file in the fastafile list into a single csv file
         #named fastaHa
         from Bio import SeqIO
         from propy import ProCheck
         import os                 
         file_to_delete = tempPath + "*.fastaHA"#txtHA"
         if os.path.exists(file_to_delete):
             os.system("rm " + file_to_delete)
         PAAC={}
         nFile = len(fastaFiles)
         for i in range(nFile):  
             f1 = open(fastaFilespath + fastaFiles[i], 'r')
             f2 = open(tempPath + fastaFiles[i] + "HA",'w')
             seqCount = self.countFastaSeq(fastaFilespath + fastaFiles[i])
             if seqCount > 1:
                  allrecord = SeqIO.parse(f1, "fasta")
                  for record in allrecord:
                      seq = str(record.seq)
                      seq = seq.upper()
                      temp = ProCheck.ProteinCheck(seq)   
                      if temp > 0:
                         if FLAG == "ENT":
                           #PAAC = fun(seq, weight=weight, maxlag=maxlag)
                           ent = self.Entropy(seq)
                           PACCstr = str(ent) 
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)
                         elif FLAG == "RENT":
                           #PAAC = fun(seq, weight=weight, maxlag=maxlag)
                           Rent = self.calRelativeEntropy(seq)
                           PACCstr = str(Rent)    
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)
                         elif FLAG == "IG":
                           IG = self.Entropy(seq) - self.calRelativeEntropy(seq)
                           PACCstr = str(IG)
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)
                         else:
                           #PAAC = fun(seq, weight=weight, maxlag=maxlag)
                           ent = self.Entropy(seq)
                           Rent = self.calRelativeEntropy(seq)
                           IG = ent - Rent
                           PACCstr = str(ent) + "," + str(Rent) + "," + str(IG)
                           #PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)


             elif seqCount == 1:
                  record = SeqIO.read(f1, "fasta")
                  seq = str(record.seq)
                  seq = seq.upper()
                  temp = ProCheck.ProteinCheck(seq)   
                  if temp > 0:
                         if FLAG == "ENT":
                           #PAAC = fun(seq, weight=weight, maxlag=maxlag)
                           ent = self.Entropy(seq)
                           PACCstr = str(ent)
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)
                         elif FLAG == "RENT":
                           #PAAC = fun(seq, weight=weight, maxlag=maxlag)
                           Rent = self.calRelativeEntropy(seq)
                           PACCstr = str(Rent)
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)
                         elif FLAG == "IG":
                           IG = self.Entropy(seq) - self.calRelativeEntropy(seq)
                           PACCstr = str(IG)
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)
                         else:
                           #PAAC = fun(seq, weight=weight, maxlag=maxlag)
                           ent = self.Entropy(seq)
                           Rent = self.calRelativeEntropy(seq)
                           IG = ent - Rent
                           PACCstr = str(ent) + "," + str(Rent) + "," + str(IG)
                           #PACCstr = str(PAAC.values()).strip('[]') + ',\n' 
                           PACCstr = PACCstr.strip('[]') + ',\n'
                           f2.write(PACCstr)

             f2.close() 
             f1.close()
             if os.path.exists(tempPath + fastaFiles[i] + "_full_features.csv") and os.path.getsize(tempPath + fastaFiles[i] + "_full_features.csv") > 0:
                
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'r+') as f4:
                        lines = f4.readlines()
                        f4.seek(0)
                        f4.truncate()
                        for line in lines:
                            
                            line = line.rstrip('\n')
                            line += f2.readline()
                            
                            f4.write(line)
             else:
                with open(tempPath + fastaFiles[i] + "HA",'r') as f2:
                    with open(tempPath + fastaFiles[i] + "_full_features.csv",'w') as f4:
                        for line in f2:
                            f4.write(line)
             f3 = open(tempPath + fastaFiles[i] + '_proteinFeaturesList.txt','a')
             des = ["Entropy", "RelEntropy","IG"]
             for PAA in des:
                   f3.write(PAA + '\n')
             f3.close()    
         #***************************************** end of entropy zone ********************



   def fetchFiles(self, mydir, ext):
       """
       Purpose: Fetch the file feature list and description 
       Argument: the directory and files extention. The file should be in one directory
       Return: File name of each class and the class description from the file name 
       """
       import os
       FileNewList=[]
       desc = []
       path = mydir
       fileList = os.listdir(path)
       for fil in fileList:
          if fil.endswith('.' + ext):
             sub = fil[0:fil.find('.')]
             desc.append(sub)
             FileNewList.append(fil)
       return FileNewList, desc 


   def convertCsvsToLabeledDataset(self, tempPath, CsvFileList, classNames, CVSFileName):
        """
        Purpose: Create SciKit dataset from features and label
        Argument: features in comma delimited and label as int
        Output: A file of sciKit dataset 
        Argument: 
             tempPath : where CsvFileList files are
             CsvFileList : the list of the single csv file created earlier for each fasta file
             CVSFileName : The name of the created single CSV file with class labels  
        """
        #delete old csv file with the same name if exists
        import os                 
        file_to_delete = CVSFileName
        if os.path.exists(file_to_delete):
            os.system("rm " + file_to_delete)

        nFile = len(CsvFileList) 
        #oFile = open('../models/Classes.lbl','w')

        FullFilePath2 = CVSFileName.split('/')
        path2=""
        for i in range(len(FullFilePath2)-1):
             path2 =path2 + str(FullFilePath2[i]) + "/"

        oFile = open(path2 + '/Classes.lbl','w')
        for i in range(nFile):
            txt = str(int(i)) + " " + (classNames[i])[0:7] + "\n" 
            oFile.write(txt)
        oFile.close()

        for i in range(nFile):
           inFileName = tempPath + "/" + CsvFileList[i]
           f1 = open(inFileName,'r')
           classID = str(int(i)) #this is the class label
           for line in f1:                      
               newLine = line.strip("\n") + classID + "\n"
               FileName = "sciKitPart" + str(int(i)) + ".csv"
               #with open(FileName,'a') as myFile: # this create a file for each class
               with open(CVSFileName,'a') as myFile:  # this create a file for all classes  
                    myFile.write(newLine)
           f1.close()
         



#if __name__ == '__main__':
    
     #pDiscrs = Proteindescriptors()
     #pProcess = PDescriptorsProcessing()

     #Hydrophobicity = pDiscrs.getHydrophobicity()
     #pK1 = pDiscrs.getpK1()
     #residuemass = pDiscrs.getResidueMass() 

     #ID = 'KRIW790103'
     #1: input fasta sequence files
     #fastaList = ['test1.fasta', 'test2.fasta','test3.fasta', 'test4.fasta']

     #2: convert each fasta sequence file to a single Csv file for a specific descriptor
     #pDiscrs.getAAComposition
     #pProcess.convertFastasToCsvkspace(pDiscrs.getGetSequenceOrderCorrelationFactorForAPAAC,'./', './', fastaList, 1)

     #3: create a single dataset (csv) with class label in the last column
     #CsvFileList, classNames = pProcess.fetchFiles('./', 'fastaHA')
     #pProcess.convertCsvsToLabeledDataset('./', CsvFileList, classNames, '/home/hamid/doc/gMine/test/pseaac/mytestFile')
     

def getFileList(dir):
   return os.listdir(dir)

def DescriptorList(Description):
       pDiscrs = Proteindescriptors()
       pPross = PDescriptorsProcessing()
       fun = pDiscrs.getSpectrumDict
       FLAG=""
       if Description =='Amino acid composition':
           fun = pDiscrs.getAAComposition
       elif Description =='Dipeptide composition':
           fun = pDiscrs.getDipeptideComp
       elif Description =='AA composition, dipeptide and 3-mers':
           fun = pDiscrs.getAADipeptideComposition
       elif Description =='Spectrum descriptors of 3-mers':
           fun = pDiscrs.getSpectrumDict

       elif Description =='Shannon entropy':
           #FLAG = "ENT"
            fun = "ENT"
       elif Description =='Relative entropy':
           #FLAG = "RENT"
           fun = "RENT"
       elif Description =='Information gain':
           #FLAG = "IG"
           fun = "IG"
       elif Description =='Entropy, relative entropy, and gain':
           #FLAG = "ERG"
           fun = "ERG"

       elif Description =='Conjoint triad':
           fun = pDiscrs.CalculateConjointTriad
       elif Description =='Predicted solvent accessibility':
           fun = pDiscrs.processASA


       elif Description =='Composition Translation Distribution of AADs':
           fun = pDiscrs.getCTD
       elif Description =='Composition descriptors of AADs':
           fun = pDiscrs.getComposition
       elif Description =='Translation descriptors of AADs':
           fun = pDiscrs.getTranslation
       elif Description =='Distribution descriptors of AADs':
           fun = pDiscrs.getDistribution
       elif Description =='Composition descriptors based on charge of AADs':
           fun = pDiscrs.getCalculateCompositionCharge
       elif Description =='Composition descriptors based on hydrophobicity of AADs':
           fun = pDiscrs.getCalculateCompositionHydrophobicity
       elif Description =='Composition descriptors based on normalized VDWV of AADs':
           fun = pDiscrs.getCalculateCompositionNormalizedVDWV
       elif Description =='Composition descriptors based on polarity of AADs':
           fun = pDiscrs.getCalculateCompositionPolarity
       elif Description =='Composition descriptors based on polarizability of AADs':
           fun = pDiscrs.getCalculateCompositionPolarizability
       elif Description =='Composition descriptors based on 2ndary structure of AADs':
           fun = pDiscrs.getCalculateCompositionSecondaryStr
       elif Description =='Composition descriptors based on solvent accessibility of AADs':
           fun = pDiscrs.getCalculateCompositionSolventAccessibility
       elif Description =='Distribution descriptors based on charge of AADs':
           fun = pDiscrs.getCalculateDistributionCharge
       elif Description =='Distribution descriptors based on hydrophobicity of AADs':
           fun = pDiscrs.getCalculateDistributionHydrophobicity
       elif Description =='Distribution descriptors based on normalized VDWV of AADs':
           fun = pDiscrs.getCalculateDistributionNormalizedVDWV
       elif Description =='Distribution descriptors based on polarity of AADs':
           fun = pDiscrs.getCalculateDistributionPolarity
       elif Description =='Distribution descriptors based on polarizability of AADs':
           fun = pDiscrs.getCalculateDistributionPolarizability
       elif Description =='Distribution descriptors based on 2ndary structure of AADs':
           fun = pDiscrs.getCalculateDistributionSecondaryStr
       elif Description =='Distribution descriptors based on solvent accessibility of AADs':
           fun = pDiscrs.getCalculateDistributionSolventAccessibility
       elif Description =='Transition descriptors based on charge of AADs':
           fun = pDiscrs.getCalculateTransitionCharge
       elif Description =='Transition descriptors based on hydrophobicity of AADs':
           fun = pDiscrs.getCalculateTransitionHydrophobicity
       elif Description =='Transition descriptors based on normalized VDWV of AADs':
           fun = pDiscrs.getCalculateTransitionNormalizedVDWV
       elif Description =='Transition descriptors based on polarity of AADs':
           fun = pDiscrs.getCalculateTransitionPolarity
       elif Description =='Transition descriptors based on polarizability of AADs':
           fun = pDiscrs.getCalculateTransitionPolarizability
       elif Description =='Transition descriptors based on 2ndary structure of AADs':
           fun = pDiscrs.getCalculateTransitionSecondaryStr
       elif Description =='Transition descriptors based on solvent accessibility of AADs':
           fun = pDiscrs.getCalculateTransitionSolventAccessibility


       elif Description =='Geary autocorrelation':
           fun = pDiscrs.getGearyAutoCorr
       elif Description =='Custom Geary autocorrelation with ID':
           fun = pDiscrs.getcusGearyAutop1  #(ID)
       elif Description =='Custom Geary autocorrelation with AAP':
           fun = pDiscrs.getcusGearyAutop2  #(AAP)
       elif Description =='Moran autocorrelation descriptors':
           fun = pDiscrs.getGetMoranAuto
       elif Description =='Moran autocorrelation descriptors with ID':
           fun = pDiscrs.getGetMoranAutop1 #(ID)
       elif Description =='Moran autocorrelation descriptors with AAP':
           fun = pDiscrs.getGetMoranAutop2 #(AAP)
       elif Description =='Normalized Moreau-Broto autocorrelation':
           fun = pDiscrs.getGetMoreauBrotoAuto
       elif Description =='Normalized Moreau-Broto autocorrelation with ID':
           fun = pDiscrs.getGetMoreauBrotoAutop1  #(ID)
       elif Description =='Normalized Moreau-Broto autocorrelation with AAP':
           fun = pDiscrs.getGetMoreauBrotoAutop2 #(AAP)

       elif Description =='Type I Pseudo amino acid composition':
           fun = pDiscrs.getGetPAAC #(lamda=10, weight=0.05)
       elif Description =='Type I Pseudo amino acid composition with ID':
           fun = pDiscrs.getGetPAACp1 #(ID, lamda=10, weight=0.05)
       elif Description =='Type I Pseudo amino acid composition with AAP':
           fun = pDiscrs.getGetPAACp2  #(AAP, lamda=10, weight=0.05)
       elif Description =='Type II pseudo-amino acid compostion':
           fun = pDiscrs.getGetAPseudoAAC1  #lamda=10, weight=0.5)

       elif Description =='Sequence order correlation factor':
           fun = pDiscrs.getGetSequenceOrderCorrelationFactor #(ID1, ID2, k=1)
       elif Description =='Sequence order correlation factor [H-phob.& H-philic.]':
           fun = pDiscrs.getGetSequenceOrderCorrelationFactorForAPAAC #(k=1) 
       elif Description =='Quasi sequence order descriptors':
           fun = pDiscrs.getGetQSO  #(maxlag=30, weight=0.1)
       elif Description =='Sequence order coupling numbers':
           fun = pDiscrs.getGetSOCN  #(maxlag=45)
       elif Description =='Quasi sequence order descriptors (matrix)':
           fun = pDiscrs.getGetQSOp #(distancematrix, maxlag=30, weight=0.1)
       elif Description =='Sequence order coupling numbers (matrix)':
           fun = pDiscrs.getGetSOCNp  #(distancematrix, maxlag=45)
       elif Description =='k-Spaced Amino Acid Pairs':
           fun = pDiscrs.getCKSAAP
       return fun


def createCSV(uploadFolder, myTempPath, FeaturesTypeRDO, ComboAAC, ComboAAD,ComboACD,\
              ComboPAAC, ComboSOC, ComboENT, ComboOTH, ComboAAP, Lambda, Weight, maxLag, ComboMatrix, aaIndex, idNum,\
              idNum2, AAproperty, kspace, email, sugFileName, chkUserAAP, AAP2):

     #uploadedFilePath = os.path.join(settings.MEDIA_ROOT, uploadFolder)
     uploadedFilePath = uploadFolder
     FileList = getFileList(uploadedFilePath)

     if len(FileList) > 0:
       #os.system("rm " + myTempPath + "/*.*")

       #the title
       if FeaturesTypeRDO =="AAC":
          title = ComboAAC
       elif FeaturesTypeRDO == "AAD":
          title = ComboAAD
       elif FeaturesTypeRDO == "ACD":
          title = ComboACD
       elif FeaturesTypeRDO == "PAAC":
          title = ComboPAAC
       elif FeaturesTypeRDO == "SOC":
          title = ComboSOC
       elif FeaturesTypeRDO == "ENT":
          title = ComboENT
       elif FeaturesTypeRDO == "OTH":
          title = ComboOTH
       elif FeaturesTypeRDO == "AAP":
          title = ComboAAP


       #csv file name
       CSVFileName = myTempPath + "/" + sugFileName
       #Step 1: Fasta file list
       path = uploadedFilePath + "/"
       #FileList =[] # see above

       #2: convert each fasta sequence file to a single Csv file for a specific descriptor
       #AAP = {}
       pDiscrs = Proteindescriptors()
       pProcess = PDescriptorsProcessing()
       if FeaturesTypeRDO=="ACD":
         if ComboACD=='Custom Geary autocorrelation with ID' \
               or ComboACD=='Moran autocorrelation descriptors with ID' \
               or ComboACD=='Normalized Moreau-Broto autocorrelation with ID':
           ID = idNum
           fun = DescriptorList(ComboACD)
           pProcess.convertFastasToCsvID(fun, path, myTempPath + '/', FileList, ID)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('ID ' + ID)
           optFile.close()
         elif ComboACD=='Custom Geary autocorrelation with AAP'\
               or ComboACD=='Moran autocorrelation descriptors with AAP' \
               or ComboACD=='Normalized Moreau-Broto autocorrelation with AAP':
           AAP = {}
           if chkUserAAP == 1:
              AAP = AAP2
              """
              AAP = {'A': A, 'C': C, 'E': E, 'D': D, 'G': G, \
                     'I': I, 'H': H, 'K': K, 'M': M, 'L': L, \
                     'N': N, 'Q': Q, 'P': P, 'S': S, 'R': R, \
                     'T': T, 'W': W, 'V': V, 'Y': Y}
              """
              fun = DescriptorList(ComboACD)
              try:
                  pProcess.convertFastasToCsvAAP(fun, path, myTempPath + "/", FileList, AAP,'')
              except ZeroDivisionError:
                  pass
           else:
             if AAproperty=='Hydrophobicity':
               AAP =  pDiscrs.getHydrophobicity()
               pname = 'Hydrophob'
             elif AAproperty=='pK1':
               AAP =  pDiscrs.getpK1()
               pname = 'pk1'
             elif AAproperty=='Residue mass':
               AAP =  pDiscrs.getResidueMass()
               pname = 'ResMass'
             elif AAproperty=='Hydrophilicity':
               AAP =  pDiscrs.getHydrophilicity()
               pname = 'Hydrophil'
             elif AAproperty=='Isoelectric point (pI)':
               AAP =  pDiscrs.getPI()
               pname = 'pI'
             elif AAproperty=='pK2':
               AAP =  pDiscrs.getpK2()
               pname = 'pk2'
             fun = DescriptorList(ComboACD)
             pProcess.convertFastasToCsvAAP(fun, path, myTempPath + "/", FileList, AAP, pname)
             #Create a file for options

           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('AAP #' + str(AAP))
           optFile.close()
         else:
           fun = DescriptorList(ComboACD)
           pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.close()



       elif FeaturesTypeRDO=="PAAC":
         if ComboPAAC=='Type I Pseudo amino acid composition' \
               or ComboPAAC=='Type II pseudo-amino acid compostion':
           lamda = int(Lambda)
           weight = float(Weight)
           fun = DescriptorList(ComboPAAC)
           pProcess.convertFastasToCsvLamdaWeight(fun, path, myTempPath + "/", FileList, lamda, weight)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('lamda ' + str(lamda) + '\n')
           optFile.write('weight ' + str(weight))
           optFile.close()
         elif ComboPAAC=='Type I Pseudo amino acid composition with ID':
           lamda = int(Lambda)
           weight = float(Weight)
           ID = str(idNum)
           fun = DescriptorList(ComboPAAC)
           pProcess.convertFastasToCsvLamdaWeightID(fun, path, myTempPath + "/", FileList, lamda, weight, ID)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('lamda ' + str(lamda) + '\n')
           optFile.write('weight ' + str(weight) + '\n')
           optFile.write('ID ' + ID)
           optFile.close()
         elif ComboPAAC=='Type I Pseudo amino acid composition with AAP':
                 #or ComboPAAC=='Type II pseudo-amino acid compostion':
           lamda = int(Lambda)
           weight = float(Weight)
           AAP = {}
           if chkUserAAP ==1:
              AAP = AAP2
              """
              AAP = {'A': A, 'C': C, 'E': E, 'D': D, 'G': G, \
                     'I': I, 'H': H, 'K': K, 'M': M, 'L': L, \
                     'N': N, 'Q': Q, 'P': P, 'S': S, 'R': R, \
                     'T': T, 'W': W, 'V': V, 'Y': Y}
              """
              fun = DescriptorList(ComboACD)
              try:
                 pProcess.convertFastasToCsvAAP(fun, path, myTempPath + "/", FileList, AAP)
              except ZeroDivisionError:
                  pass

           else:
             if AAproperty=='Hydrophobicity':
               AAP =  [pDiscrs.getHydrophobicity()]
               pname = 'Hydrophob'
             elif AAproperty=='pK1':
               AAP =  [pDiscrs.getpK1()]
               pname = 'pk1'
             elif AAproperty=='Residue mass':
               AAP =  [pDiscrs.getResidueMass()]
               pname = 'ResMass'
             elif AAproperty=='Hydrophilicity':
               AAP =  [pDiscrs.getHydrophilicity()]
               pname = 'Hydrophil'
             elif AAproperty=='Isoelectric point (pI)':
               AAP =  [pDiscrs.getPI()]
               pname = 'pI'
             elif AAproperty=='pK2':
               AAP =  [pDiscrs.getpK2()]
               pname = 'pk2'
            
             fun = DescriptorList(ComboPAAC)
             pProcess.convertFastasToCsvLamdaWeightAAP(fun, path, myTempPath + "/", FileList, lamda, weight, AAP, pname)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('lamda ' + str(lamda) + '\n')
           optFile.write('weight ' + str(weight) + '\n')
           optFile.write('AAP #' + str(AAP))
           optFile.close()

       elif FeaturesTypeRDO=="SOC":
         if ComboSOC=='Quasi sequence order descriptors':
           weight = float(Weight)
           maxlag = int(maxLag)
           fun = DescriptorList(ComboSOC)
           pProcess.convertFastasToCsvWeightMaxlag(fun, path, myTempPath + "/", FileList, weight, maxlag)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('maxlag ' + str(maxlag) + '\n')
           optFile.write('weight ' + str(weight))
           optFile.close()
         elif ComboSOC=='Quasi sequence order descriptors (matrix)':
           weight = float(Weight)
           maxlag = int(maxLag)
           mat = ''
           if ComboMatrix =="Grantham distance matrix":
              matrixFile = "mats/Grantham.csv"
              mat = 'Gran'
           else:
              matrixFile = "mats/Schneider-Wrede.csv"
              mat = 'SW'
           matrix = pDiscrs.distMatrix(matrixFile)
           fun = DescriptorList(ComboSOC)
           pProcess.convertFastasToCsvWeightMaxlagMatrix(fun, path, myTempPath + "/", FileList, weight, maxlag, matrix, mat)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('maxlag ' + str(maxlag) + '\n')
           optFile.write('weight ' + str(weight) + '\n')
           optFile.write('matrix ' + ComboMatrix + "\n")
           optFile.close()

         elif ComboSOC=='Sequence order coupling numbers':
           maxlag = int(maxLag)
           fun = DescriptorList(ComboSOC)
           pProcess.convertFastasToCsvMaxlag(fun, path, myTempPath + "/", FileList, maxlag)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('maxlag ' + str(maxlag))
           optFile.close()

         elif ComboSOC=='Sequence order coupling numbers (matrix)':
           maxlag = int(maxLag)
           if ComboMatrix =="Grantham distance matrix":
              matrixFile = "mats/Grantham.csv"
              mat = 'Gran'
           else:
              matrixFile = "mats/Schneider-Wrede.csv"
              mat = 'SW'
           matrix = pDiscrs.distMatrix(matrixFile)
           fun = DescriptorList(ComboSOC)
           pProcess.convertFastasToCsvMaxlagMatrix(fun, path, myTempPath + "/", FileList, maxlag, matrix, mat)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('maxlag ' + str(maxlag) + '\n')
           optFile.write('matrix ' + ComboMatrix)
           optFile.close()




         elif ComboSOC=='Sequence order correlation factor':
           if Lambda is not None:
               lamda = int(Lambda)
           ID1 = str(idNum)
           ID2 = str(idNum2)
           fun = DescriptorList(ComboSOC)
           pProcess.convertFastasToCsv2IDs(fun, path, myTempPath + '/', FileList, ID1, ID2)
           #Create a file for options
           optFile = open(myTempPath + "/" + 'options.txt','a')
           optFile.write('lamda ' + str(lamda) + '\n')
           optFile.write('ID1 #' + ID1)
           optFile.write('ID1 #' + ID2)
           optFile.close()
         elif ComboSOC=='Sequence order correlation factor [H-phob.& H-philic.]':
          try:
             kspace = int(kspace)
             fun = DescriptorList(ComboSOC)
             pProcess.convertFastasToCsvkspace(fun, path, myTempPath + '/', FileList, kspace)
             #Create a file for options
             optFile = open(myTempPath + "/" + 'options.txt','a')
             optFile.write('kspace ' + str(kspace))
             optFile.close()
          except KeyError:
             pass
         else:
           if ComboACD:
             fun = DescriptorList(ComboACD)
             pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
             #Create a file for options
             optFile = open(myTempPath + "/" + 'options.txt','a')
             optFile.write('None')
             optFile.close()
             """
             elif ComboPAAC:
             fun = DescriptorList(ComboPAAC)
             pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
             #Create a file for options
             optFile = open(myTempPath + "/" + 'options.txt','w')
             optFile.write('None')
             optFile.close()
             """
           elif ComboSOC:
             fun = DescriptorList(ComboSOC)
             pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
             #Create a file for options
             optFile = open(myTempPath + "/" + 'options.txt','a')
             optFile.write('None')
             optFile.close()


       elif FeaturesTypeRDO=='AAC':
            fun = DescriptorList(ComboAAC)
            pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
            #Create a file for options
            optFile = open(myTempPath + "/" + 'options.txt','a')
            optFile.write('None')
            optFile.close()
            
       elif FeaturesTypeRDO=='AAP':
            fun = DescriptorList(ComboAAP)
            pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
            #Create a file for options
            optFile = open(myTempPath + "/" + 'options.txt','a')
            optFile.write('None')
            optFile.close()

       elif FeaturesTypeRDO=='AAD':
            fun = DescriptorList(ComboAAD)
            pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
            #Create a file for options
            optFile = open(myTempPath + "/" + 'options.txt','a')
            optFile.write('None')
            optFile.close()

       elif FeaturesTypeRDO=='ENT':
            fun = DescriptorList(ComboENT)
            #if fun == "ENT"
            pProcess.convertFastasToCsvEntropy_relEntropy(path, myTempPath + '/', FileList, fun)
            #Create a file for options
            optFile = open(myTempPath + "/" + 'options.txt','a')
            optFile.write('None')
            optFile.close()
       elif FeaturesTypeRDO=='OTH':
            if ComboOTH=='Conjoint triad':
               fun = DescriptorList(ComboOTH)
               pProcess.convertFastasToCsv0(fun, path, myTempPath + '/', FileList)
               #Create a file for options
               optFile = open(myTempPath + "/" + 'options.txt','a')
               optFile.write('None')
               optFile.close()
            else:
               fun = DescriptorList(ComboOTH)
               pProcess.convertFastasToCsvASA(fun, path, myTempPath + '/', FileList)
               #Create a file for options
               optFile = open(myTempPath + "/" + 'options.txt','a')
               optFile.write('None')
               optFile.close()



     else:
       pass 
def seqCorr(idN):
    id1,id2 = idN
    createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = 'SOC', ComboAAC = '', ComboAAD = '',\
        ComboACD = '', ComboPAAC = '', ComboSOC = 'Sequence order correlation factor', ComboENT = '', ComboOTH = '', ComboAAP = '',\
        Lambda = Lambda, Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',\
        idNum = id1,idNum2 = id2, AAproperty = '', kspace = '',  email = '',\
        sugFileName = 'S_Features.csv', chkUserAAP = '', AAP2 = '')

def feps(infile, outfile):#main(argv):
    global uploadFolder
    global myTempPath
    global Lambda
    uploadFolder = infile#''
    myTempPath = outfile#''


    problem_features = []
    FeaturesTypeRDO = ['AAP','OTH','ENT','AAD','AAC','ACD','PAAC','SOC']
    ComboAAC = ['Dipeptide composition','Spectrum descriptors of 3-mers']
    ComboAAD = ['Composition Translation Distribution of AADs','Composition descriptors of AADs',\
                'Translation descriptors of AADs','Distribution descriptors of AADs',\
                'Composition descriptors based on charge of AADs',\
                'Composition descriptors based on hydrophobicity of AADs',\
                'Composition descriptors based on normalized VDWV of AADs',\
                'Composition descriptors based on polarity of AADs',\
                'Composition descriptors based on polarizability of AADs',\
                'Composition descriptors based on 2ndary structure of AADs',\
                'Composition descriptors based on solvent accessibility of AADs',\
                'Distribution descriptors based on charge of AADs',\
                'Distribution descriptors based on hydrophobicity of AADs',\
                'Distribution descriptors based on normalized VDWV of AADs',\
                'Distribution descriptors based on polarity of AADs',\
                'Distribution descriptors based on polarizability of AADs',\
                'Distribution descriptors based on 2ndary structure of AADs',\
                'Distribution descriptors based on solvent accessibility of AADs',\
                'Transition descriptors based on charge of AADs','Transition descriptors based on hydrophobicity of AADs',\
                'Transition descriptors based on normalized VDWV of AADs','Transition descriptors based on polarity of AADs',\
                'Transition descriptors based on polarizability of AADs',\
                'Transition descriptors based on 2ndary structure of AADs',\
                'Transition descriptors based on solvent accessibility of AADs']

    ComboACD = ['Geary autocorrelation','Moran autocorrelation descriptors','Normalized Moreau-Broto autocorrelation']

    ComboPAAC = ['Type I Pseudo amino acid composition','Type II pseudo-amino acid compostion']
    ComboSOC = ['Quasi sequence order descriptors','Quasi sequence order descriptors (matrix)',\
                'Sequence order coupling numbers','Sequence order coupling numbers (matrix)',\
                'Sequence order correlation factor [H-phob.& H-philic.]']    
    ComboENT = 'Entropy, relative entropy, and gain'
    ComboOTH = 'Conjoint triad'
    ComboAAP = 'k-Spaced Amino Acid Pairs'
    Lambda = 30 
    Weight = 0.05 #0.1 For SOC
    MaxLag = 30
    ComboMatrix = ['Grantham distance matrix','Schneider-Wrede distance matrix']      
    idNum = ['ANDN920101','ARGP820101','ARGP820102','ARGP820103','BEGF750101','BEGF750102','BEGF750103','BHAR880101',\
             'BIGC670101','BIOV880101','BIOV880102','BROC820101','BROC820102','BULH740101','BULH740102','BUNA790101',\
             'BUNA790102','BUNA790103','BURA740101','BURA740102','CHAM810101','CHAM820101','CHAM820102','CHAM830101',\
             'CHAM830102','CHAM830103','CHAM830104','CHAM830105','CHAM830106','CHAM830107','CHAM830108','CHOC750101',\
             'CHOC760101','CHOC760102','CHOC760103','CHOC760104','CHOP780101','CHOP780201','CHOP780202','CHOP780203',\
             'CHOP780204','CHOP780205','CHOP780206','CHOP780207','CHOP780208','CHOP780209','CHOP780210','CHOP780211',\
             'CHOP780212','CHOP780213','CHOP780214','CHOP780215','CHOP780216','CIDH920101','CIDH920102','CIDH920103',\
             'CIDH920104','CIDH920105','COHE430101','CRAJ730101','CRAJ730102','CRAJ730103','DAWD720101','DAYM780101',\
             'DAYM780201','DESM900101','DESM900102','EISD840101','EISD860101','EISD860102','EISD860103','FASG760101',\
             'FASG760102','FASG760103','FASG760104','FASG760105','FAUJ830101','FAUJ880101','FAUJ880102','FAUJ880103',\
             'FAUJ880104','FAUJ880105','FAUJ880106','FAUJ880107','FAUJ880108','FAUJ880109','FAUJ880110','FAUJ880111',\
             'FAUJ880112','FAUJ880113','FINA770101','FINA910101','FINA910102','FINA910103','FINA910104','GARJ730101',\
             'GEIM800101','GEIM800102','GEIM800103','GEIM800104','GEIM800105','GEIM800106','GEIM800107','GEIM800108',\
             'GEIM800109','GEIM800110','GEIM800111','GOLD730101','GOLD730102','GRAR740101','GRAR740102','GRAR740103',\
             'GUYH850101','HOPA770101','HOPT810101','HUTJ700101','HUTJ700102','HUTJ700103','ISOY800101','ISOY800102',\
             'ISOY800103','ISOY800104','ISOY800105','ISOY800106','ISOY800107','ISOY800108','JANJ780101','JANJ780102',\
             'JANJ780103','JANJ790101','JANJ790102','JOND750101','JOND750102','JOND920101','JOND920102','JUKT750101',\
             'JUNJ780101','KANM800101','KANM800102','KANM800103','KANM800104','KARP850101','KARP850102','KARP850103',\
             'KHAG800101','KLEP840101','KRIW710101','KRIW790101','KRIW790102','KRIW790103','KYTJ820101','LAWE840101',\
             'LEVM760101','LEVM760102','LEVM760103','LEVM760104','LEVM760105','LEVM760106','LEVM760107','LEVM780101',\
             'LEVM780102','LEVM780103','LEVM780104','LEVM780105','LEVM780106','LEWP710101','LIFS790101','LIFS790102',\
             'LIFS790103','MANP780101','MAXF760101','MAXF760102','MAXF760103','MAXF760104','MAXF760105','MAXF760106',\
             'MCMT640101','MEEJ800101','MEEJ800102','MEEJ810101','MEEJ810102','MEIH800101','MEIH800102','MEIH800103',\
             'MIYS850101','NAGK730101','NAGK730102','NAGK730103','NAKH900101','NAKH900102','NAKH900103','NAKH900104',\
             'NAKH900105','NAKH900106','NAKH900107','NAKH900108','NAKH900109','NAKH900110','NAKH900111','NAKH900112',\
             'NAKH900113','NAKH920101','NAKH920102','NAKH920103','NAKH920104','NAKH920105','NAKH920106','NAKH920107',\
             'NAKH920108','NISK800101','NISK860101','NOZY710101','OOBM770101','OOBM770102','OOBM770103','OOBM770104',\
             'OOBM770105','OOBM850101','OOBM850102','OOBM850103','OOBM850104','OOBM850105','PALJ810101','PALJ810102',\
             'PALJ810103','PALJ810104','PALJ810105','PALJ810106','PALJ810107','PALJ810108','PALJ810109','PALJ810110',\
             'PALJ810111','PALJ810112','PALJ810113','PALJ810114','PALJ810115','PALJ810116','PARJ860101','PLIV810101',\
             'PONP800101','PONP800102','PONP800103','PONP800104','PONP800105','PONP800106','PONP800107','PONP800108',\
             'PRAM820101','PRAM820102','PRAM820103','PRAM900101','PRAM900102','PRAM900103','PRAM900104','PTIO830101',\
             'PTIO830102','QIAN880101','QIAN880102','QIAN880103','QIAN880104','QIAN880105','QIAN880106','QIAN880107',\
             'QIAN880108','QIAN880109','QIAN880110','QIAN880111','QIAN880112','QIAN880113','QIAN880114','QIAN880115',\
             'QIAN880116','QIAN880117','QIAN880118','QIAN880119','QIAN880120','QIAN880121','QIAN880122','QIAN880123',\
             'QIAN880124','QIAN880125','QIAN880126','QIAN880127','QIAN880128','QIAN880129','QIAN880130','QIAN880131',\
             'QIAN880132','QIAN880133','QIAN880134','QIAN880135','QIAN880136','QIAN880137','QIAN880138','QIAN880139',\
             'RACS770101','RACS770102','RACS770103','RACS820101','RACS820102','RACS820103','RACS820104','RACS820105',\
             'RACS820106','RACS820107','RACS820108','RACS820109','RACS820110','RACS820111','RACS820112','RACS820113',\
             'RACS820114','RADA880101','RADA880102','RADA880103','RADA880104','RADA880105','RADA880106','RADA880107',\
             'RADA880108','RICJ880101','RICJ880102','RICJ880103','RICJ880104','RICJ880105','RICJ880106','RICJ880107',\
             'RICJ880108','RICJ880109','RICJ880110','RICJ880111','RICJ880112','RICJ880113','RICJ880114','RICJ880115',\
             'RICJ880116','RICJ880117','ROBB760101','ROBB760102','ROBB760103','ROBB760104','ROBB760105','ROBB760106',\
             'ROBB760107','ROBB760108','ROBB760109','ROBB760110','ROBB760111','ROBB760112','ROBB760113','ROBB790101',\
             'ROSG850101','ROSG850102','ROSM880101','ROSM880102','ROSM880103','SIMZ760101','SNEP660101','SNEP660102',\
             'SNEP660103','SNEP660104','SUEM840101','SUEM840102','SWER830101','TANS770101','TANS770102','TANS770103',\
             'TANS770104','TANS770105','TANS770106','TANS770107','TANS770108','TANS770109','TANS770110','VASM830101',\
             'VASM830102','VASM830103','VELV850101','VENT840101','VHEG790101','WARP780101','WEBA780101','WERD780101',\
             'WERD780102','WERD780103','WERD780104','WOEC730101','WOLR810101','WOLS870101','WOLS870102','WOLS870103',\
             'YUTK870101','YUTK870102','YUTK870103','YUTK870104','ZASB820101','ZIMJ680101','ZIMJ680102','ZIMJ680103',\
             'ZIMJ680104','ZIMJ680105','AURR980101','AURR980102','AURR980103','AURR980104','AURR980105','AURR980106',\
             'AURR980107','AURR980108','AURR980109','AURR980110','AURR980111','AURR980112','AURR980113','AURR980114',\
             'AURR980115','AURR980116','AURR980117','AURR980118','AURR980119','AURR980120','ONEK900101','ONEK900102',\
             'VINM940101','VINM940102','VINM940103','VINM940104','MUNV940101','MUNV940102','MUNV940103','MUNV940104',\
             'MUNV940105','WIMW960101','KIMC930101','MONM990101','BLAM930101','PARS000101','PARS000102','KUMS000101',\
             'KUMS000102','KUMS000103','KUMS000104','TAKK010101','FODM020101','NADH010101','NADH010102','NADH010103',\
             'NADH010104','NADH010105','NADH010106','NADH010107','MONM990201','KOEP990101','KOEP990102','CEDJ970101',\
             'CEDJ970102','CEDJ970103','CEDJ970104','CEDJ970105','FUKS010101','FUKS010102','FUKS010103','FUKS010104',\
             'FUKS010105','FUKS010106','FUKS010107','FUKS010108','FUKS010109','FUKS010110','FUKS010111','FUKS010112',\
             'AVBF000101','AVBF000102','AVBF000103','AVBF000104','AVBF000105','AVBF000106','AVBF000107','AVBF000108',\
             'AVBF000109','YANJ020101','MITS020101','TSAJ990101','TSAJ990102','COSI940101','PONP930101','WILM950101',\
             'WILM950102','WILM950103','WILM950104','KUHL950101','GUOD860101','JURD980101','BASU050101','BASU050102',\
             'BASU050103','SUYM030101','PUNT030101','PUNT030102','GEOR030101','GEOR030102','GEOR030103','GEOR030104',\
             'GEOR030105','GEOR030106','GEOR030107','GEOR030108','GEOR030109','ZHOH040101','ZHOH040102','ZHOH040103',\
             'BAEK050101','HARY940101','PONJ960101','DIGM050101','WOLR790101','OLSK800101','KIDA850101','GUYH850102',\
             'GUYH850103','GUYH850104','GUYH850105','ROSM880104','ROSM880105','JACR890101','COWR900101','BLAS910101',\
             'CASG920101','CORJ870101','CORJ870102','CORJ870103','CORJ870104','CORJ870105','CORJ870106','CORJ870107',\
             'CORJ870108','MIYS990101','MIYS990102','MIYS990103','MIYS990104','MIYS990105','ENGD860101','FASG890101']
    idNum2 = idNum
    AAproperty = ['Hydrophobicity','pK1','Residue mass','Hydrophilicity','Isoelectric point (pI)','pK2']
    kspace = 1 
    chkUserAAP = [0,1]
    sugFileName = 'S_Features.csv'
    for rdo in FeaturesTypeRDO:
        print(rdo)
        if rdo == 'ACD':
            for acd in ComboACD:
                if acd=='Custom Geary autocorrelation with ID' \
                    or acd=='Moran autocorrelation descriptors with ID' \
                    or acd=='Normalized Moreau-Broto autocorrelation with ID':
                        
                        for id in idNum:
                            
                            try:
                                createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                          ComboACD = acd, ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                          Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = id,\
                                          idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                                          chkUserAAP = '', AAP2 = '')
                            except:
                                problem_features.append(acd+'_'+id)
                                pass
                elif acd=='Custom Geary autocorrelation with AAP'\
                    or acd=='Moran autocorrelation descriptors with AAP' \
                    or acd=='Normalized Moreau-Broto autocorrelation with AAP':
                        
                        for aaprop in AAproperty:
                            
                            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                      ComboACD = acd, ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                      Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                                      idNum2 = '', AAproperty = aaprop, kspace = '',  email = '', sugFileName = sugFileName,\
                                      chkUserAAP = '', AAP2 = '')
                elif acd == 'Geary autocorrelation' or acd == 'Moran autocorrelation descriptors'\
                    or acd == 'Normalized Moreau-Broto autocorrelation':
                            
                            
                            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                      ComboACD = acd, ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                      Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                                      idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                                      chkUserAAP = '', AAP2 = '')                        
        elif rdo == 'PAAC':
            for paac in ComboPAAC:
                
                if paac=='Type I Pseudo amino acid composition' \
                    or paac=='Type II pseudo-amino acid compostion':
                        
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = paac, ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = Lambda, Weight = Weight, maxLag = '', ComboMatrix = '', aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = '', kspace = '',  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')   
                elif paac=='Type I Pseudo amino acid composition with ID':
                    for id in idNum:
                        try:
                            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                      ComboACD = '', ComboPAAC = paac, ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                      Lambda = Lambda, Weight = Weight, maxLag = '', ComboMatrix = '', aaIndex = '',\
                                      idNum = id,idNum2 = '', AAproperty = '', kspace = '',  email = '',\
                                      sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')
                        except:
                            
                            problem_features.append(paac+'_'+id)
                            pass
                elif paac=='Type I Pseudo amino acid composition with AAP':
                    for aaprop in AAproperty:
                            
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = paac, ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = Lambda, Weight = Weight, maxLag = '', ComboMatrix = '', aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = aaprop, kspace = '',  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')

        elif rdo == 'SOC':
            Weight = 0.1
            for soc in ComboSOC:
            
                if soc=='Quasi sequence order descriptors':
                        
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = '', ComboSOC = soc, ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = '', Weight = Weight, maxLag = MaxLag, ComboMatrix = '', aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = '', kspace = '',  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '') 
                elif soc=='Quasi sequence order descriptors (matrix)':
                    
                    for mat in ComboMatrix:
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = '', ComboSOC = soc, ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = '', Weight = Weight, maxLag = MaxLag, ComboMatrix = mat, aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = '', kspace = '',  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')
                elif soc=='Sequence order coupling numbers':
                        
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = '', ComboSOC = soc, ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = '', Weight = '', maxLag = MaxLag, ComboMatrix = mat, aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = '', kspace = '',  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')
                elif soc=='Sequence order coupling numbers (matrix)':
                    
                    for mat in ComboMatrix:
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = '', ComboSOC = soc, ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = '', Weight = '', maxLag = MaxLag, ComboMatrix = mat, aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = '', kspace = '',  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')
                elif soc=='Sequence order correlation factor':
                    
                    
                    for id,id2 in combinations(idNum,2):
                 
                                try:
                                    
                                    createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                        ComboACD = '', ComboPAAC = '', ComboSOC = soc, ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                        Lambda = Lambda, Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',\
                                        idNum = id,idNum2 = id2, AAproperty = '', kspace = '',  email = '',\
                                        sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')
                                except:
                                    
                                    problem_features.append(soc+'_'+id+'_'+id2)
                                    pass
                elif soc=='Sequence order correlation factor [H-phob.& H-philic.]':
                        createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                                  ComboACD = '', ComboPAAC = '', ComboSOC = soc, ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                  Lambda = '', Weight = '', maxLag = '', ComboMatrix = mat, aaIndex = '',\
                                  idNum = '',idNum2 = '', AAproperty = '', kspace = kspace,  email = '',\
                                  sugFileName = sugFileName, chkUserAAP = '', AAP2 = '')
        elif rdo == 'AAC':
            for aac in ComboAAC:                    
                            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = aac, ComboAAD = '',\
                                      ComboACD = '', ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                      Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                                      idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                                      chkUserAAP = '', AAP2 = '')
        elif rdo == 'AAD':
            for aad in ComboAAD:                    
                            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = aad,\
                                      ComboACD = '', ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = '',\
                                      Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                                      idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                                      chkUserAAP = '', AAP2 = '')
        elif rdo == 'ENT':
            ent = ComboENT
            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                        ComboACD = '', ComboPAAC = '', ComboSOC = '', ComboENT = ent, ComboOTH = '', ComboAAP = '',\
                        Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                        idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                        chkUserAAP = '', AAP2 = '')
        elif rdo == 'OTH':
            oth = ComboOTH
            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                        ComboACD = '', ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = oth, ComboAAP = '',\
                        Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                        idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                        chkUserAAP = '', AAP2 = '')
        elif rdo == 'AAP':
            aap = ComboAAP                    
            createCSV(uploadFolder, myTempPath, FeaturesTypeRDO = rdo, ComboAAC = '', ComboAAD = '',\
                        ComboACD = '', ComboPAAC = '', ComboSOC = '', ComboENT = '', ComboOTH = '', ComboAAP = aap,\
                        Lambda = '', Weight = '', maxLag = '', ComboMatrix = '', aaIndex = '',idNum = '',\
                        idNum2 = '', AAproperty = '', kspace = '',  email = '', sugFileName = sugFileName,\
                        chkUserAAP = '', AAP2 = '')
    f3 = open(myTempPath +'problem_Features.txt','w')
    for prob in problem_features:
        f3.write(prob + '\n')
    f3.close()  
