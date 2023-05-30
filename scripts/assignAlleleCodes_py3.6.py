#!/usr/bin/python3.6

"""
Allele Code Assignment Algorithm:	Using input config file (--config) to denote names of core loci in a cgMLST scheme, assess nearest neighbor
									distances of input allele profiles (--alleles) against legacy profiles contained in input data
									diretory (--datadir) in a hierarchical method at organism-specific thresholds defined by input prefix (--prefix)
									using tree file in data directory as a guide.
Required arguments:
	-a, --alleles:	new allele profiles (Key/StrainID in column 1, locus names as headers with allele numbers in rows beneath, starting in column 2)
	-c, --config:	text file containing names of core loci on separate lines with no header (i.e. all lines are considered locu names)
	-d, --datadir:	data directory (string:  folder where results will be saved it --no-save flag not given)
	-p, --prefix:	(string:  what will be appended in front of the Allele Code and subfolders in data directory where data is stored)

Optional args:
	--nosave:  if provided, tree and allele calls file(s) will not be saved, only results printed to terminal
	--verbose:	if provided, also print what is written to log file to terminal too
	-o, --output:  if provided, print results to file rather than terminal. Delimiter determined by extension (',' for csv, '\t' for tsv)

Version:
	2.1
	
Author:
	Grant Williams, MS (gmwilliams1@cdc.gov)
"""

# Stdlib imports
import os
import sys
import gzip
import json
import shutil
import traceback
import logging
import logging.handlers
import argparse
from functools import partial
from datetime import datetime, date
from collections import defaultdict, namedtuple

#============ GLOBAL VARIABLES =============#
prefix = ''			# organism-specific abbreviation prepended to Allele Codes and data directories
version = '2.1'		# current version of implementd algorithm
configPath = ''		# path to core locus names config file
coreLoci = []		# list of core locus names loaded from config file
newAllelesPath = ''	# path to cgMLST allele profiles (tsv or csv); set at run time
delim = ''			# delimiter to use in splitting new allele calls file from csv (",") or tsv ("\t")
outputPath = ''		# path to output file (tsv or csv --> same as input allele profile delimiter if not specified)

# maximum pairwise distance thresholds, weighted to percent shared loci
thresholds = {'CAMP':[100.0*i/1343 for i in [84, 61, 24, 14, 5, 1]],	# [6.95%, 4.54%, 1.79%, 1.04%, 0.372%, 0.0745%]
				'LMO':[100.0*i/1748 for i in [71, 51, 36, 19, 7, 1]],	# [4.06%, 2.92%, 2.06%, 1.09%, 0.400%, 0.0572%]
				'SALM':[100.0*i/3002 for i in [80, 28, 15, 7, 4, 1]],	# [2.66%, 0.933%, 0.5%, 0.233%, 0.133%, 0.0333%]
				'EC':[100.0*i/2513 for i in [77, 51, 16, 6, 1]]}		# [3.06%, 2.03%, 0.637%, 0.239%, 0.0398%]
defaultThresholds = [200, 150, 100, 50, 25, 1]
minpres = 0.95	# minimum percent core called (i.e. total loci with alleles > 0)
xCodeList = []	# list of Allele Codes whose within-code distance exceeds 4x the corresponding threshold
nosave = False	# whether or not to cancel overwriting files in data directory after processing
verbose = False	# whether messages sent to log file will also be sent to terminal

#============ DATA DIRECTORIES ============#
# primary data directory, set as input argument
DATA_DIR = ''

# subdirectories to be created or validated for existence
DATA_DIRS = ['tree', 
			 os.path.join('tree', 'current'), 
			 'allele_calls', 
			 os.path.join('allele_calls', 'current')]


#============ LOG DIRECTORIES ============#
# primary log directory, holding other subdirectories within
LOG_DIR = ''

# subdirectories containing other useful things
LOG_DIRS = ['change_log', 'Xcodes']


#========= OPERATIONAL VARIABLES =========#
# Counter for how many distances have been calculated
cntDistancesCalculated = 0

# Variables to track changed codes
numChanged = 0
changedKeys = {}

#========== LOGGING FUNCTIONALITY ==========#
# global logger:
Logger = None
def initialize_logging(log_dir):
	"""
	initialize_logging:  create directories for log files if not already done
		Arguments:
			log_dir:  string --> primary directory into which _nomenclature_logs folder will be made
	"""
	file_name = os.path.join(log_dir, 
							'{}_nomenclature_logs'.format(prefix),
							'wgst_log_{}.txt'.format(datetime.now().strftime("%Y-%m-%d@%H-%M-%S")))
	
	if not os.path.exists(os.path.dirname(file_name)):
		os.makedirs(os.path.dirname(file_name))
	
	global Logger

	logging.basicConfig(format='%(asctime)s\t%(message)s',
						datefmt='%m-%d-%Y %H:%M:%S',
						filename=file_name,
						level=logging.DEBUG)

	Logger = logging.getLogger()

def make_pretty(count):
	"""
	make_pretty:  return input number of tabs as a spacer
		Arguments:
			count:  int --> total number of spaces to return
		Returns:
			string
	"""
	return '\t'*count

def flush_handlers(handlers):
	"""
	flush_handlers:  clear buffers for Logger class to wipe any stored text
		Arguments:
			handlers:  list --> list of handlers to flush
	"""
	for h in handlers:
		h.flush()

def close_handlers(handlers):
	"""
	close_handlers:  close all open Logger handlers in input list
		Arguments:
			handlers:  list --> list of open Logger handlers to close
	"""
	for h in handlers:
		h.close()

def log(log_type, message, depth):
	"""
	log:  adds logging functions to Logger class if not already there, or calls them otherwise
		Arguments:
			log_type:  string --> type of message to log (info, error, or exception)
			message:   string --> text to write to log file
			depth:	   int --> total number of tabs to add to front of message text
	"""
	global Logger

	if hasattr(Logger, log_type):
		getattr(Logger, log_type)(make_pretty(depth)+message)
	else:
		getattr(Logger, 'info')(make_pretty(depth)+message)
	
	# clear logging handler buffer text
	flush_handlers(Logger.handlers)

def log_message(message, depth=0):
	"""
	log_message:  write simple message to log file.  Also print to screen if --verbose flag entered by user
		Arguments:
			message:  string --> text to write to log file
			depth:	  int --> number of tabs to add to front of message text
	"""
	log('info', message, depth)
	# also print to terminal if user opted for verbose processing
	if verbose:
		print('{}{}'.format('\t'*depth, message))
	
def log_error(message, depth=0):
	"""
	log_error:  write input error message to log file
		Arguments:
			message:  string --> text to write to log file
			depth:	  int --> number of tabs to add to front of message text
	"""
	log('error', message, depth)
	# also print to terminal if user opted for verbose processing
	if verbose:
		print('{}{}'.format('\t'*depth, message))

def log_exception(message, depth=0):
	"""
	log_exception:  write input exception message to log file
		Arguments:
			message:  string --> text to write to log file
			depth:	  int --> number of tabs to add to front of message text
	"""
	log('exception', message, depth)
	# also print to terminal if user opted for verbose processing
	if verbose:
		print('{}{}'.format('\t'*depth, message))

	
#========================= SOME TOOLING ======================================#
def initialize_srcfiles(data_dir):
	"""
	initialize_srcfiles:  make *_nomenclature_srcfiles in input data directory if not already present
		Arguments:
			data_dir:  string --> path to folder where data directory should be created
	"""
	path = os.path.join(DATA_DIR, '{}_nomenclature_srcfiles'.format(prefix))
	if not os.path.exists(path):
		os.makedirs(path)
	
def Now():
	"""
	Now:  return current time in yyyy-mm-dd hh:min:ss format
	"""
	return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def LoadProfilesFromFile(path):
	"""
	LoadProfilesFromFile: returns dictionary of profiles from input file (csv or tsv) in correct naming order
							(order pre-determined in coreLoci global variable)
		Arguments:
			path:  string --> file path to allele profiles
		Returns:
			dict --> key:[alleles] for each key (column 1) and profiles (column 2+)
	"""
	profiles = {}
	fields = []
	# get raw values
	with open(path, 'r') as f:
		# load all profiles and headers into list of lists
		content = [line.strip().split(delim) for line in f if len(line)]
		
		# get field headings
		fields = content[0][1:]
		
		# convert profiles from lists to locus name:allele dict in profile dict
		profiles = {content[i+1][0]:{fields[f]:int(content[i+1][f+1]) for f in range(len(fields))} for i in range(len(content[1:]))}
		
		# convert profiles to lists in order of global coreloci list, filling in missing loci with allele 0
		for k, (key, values) in enumerate(profiles.items()):
			profiles[key] = [values.get(locus, 0) for locus in coreLoci]
			
	# return resulting dict
	return profiles
	

#====================== Xcode functionality ===========================#
def SetXcodeList(xCodeListPath):
	"""
	SetXcodeList: set global xCodeList variable using Xcodes.tsv file in log directory if present
		Arguments:
			xCodeListPath:  string --> path to Xcodes.tsv file
	"""
	global xCodeList
	if not os.path.exists(xCodeListPath):
		# set xCodeList to empty list if no Xcodes.tsv file found
		log_message("Xcodes.tsv file not found", depth=0)
		xCodeList = []
		return
	
	try:
		with open(xCodeListPath, 'r') as f:
			# otherwise, read file into list, and update global xCodeList variable
			xFileLines = [line for line in f if len(line)]
			f.close()
			log_message("{} lines in Xcodes file".format(len(xFileLines)), depth=0)
			
			# get rid of any field headers if present
			if '.' not in xFileLines[0]:
				xFileLines.pop(0)
			
			# set xCodeList to first index of tab-delimited lines of Xcode file
			xCodeList = [line.split('\t')[0] for line in xFileLines]
			
	except:
		# Log if an error happened loading Xcodes and set xCodeList to an empty list
		log_message("Error processing Xcodes.tsv file.", depth=1)
		xCodeList = []
		
def CheckXcodeList(name):
	"""
	CheckXcodeList:  if input Allele Code matches a current Xcode up to one of its digits, then return that Xcode; return input Allele Code otherwise
		Arguments:
			name:  string --> newly assigned Allele Code for assessment on X status
		Returns:
			string (either matching Xcode or input name if no matches found)
	"""
	# this is the length of what comes in front of the real name, like: 'SALM1.0 - '
	l = len(prefix) + len(version) + 3
	# trim that off and convert input name to a list
	nameAsList = list(map(str, name[l:].split('.')))
		
	for xCode in xCodeList:
		# skip if entry's name is shorter
		if len(nameAsList) < len(xCode.split('.')):
			continue
		# if the name matches Xcode name, add an x to the end of it at the teriminal position of the Xcode
		if '.'.join(nameAsList[:len(xCode.split('.'))]) == xCode:
			return '{}{} - {}x'.format(prefix, version, xCode)
	
	# return input name otherwise
	return name
	
#====================== Class definitions =============================#
Cluster = namedtuple('Cluster', ['id', 'members', 'diameter', 'preferred'])

class Tree(object):
	"""
	Tree:  class that's essentially a dictionary with algorithm-specific member variables and functions
	"""
	CURRENT = None
	DEPTH = -1

	def __init__(self, depth):
		"""
		Initialize with default values
		"""
		Tree.CURRENT = self					# self reference
		Tree.DEPTH = depth					# total number of digits in full-length Allele Code
		self._tree = Node(1, 0, None)		# initial empty node of any Tree
		self._names = defaultdict(list)		# full-length Allele Code that will hold the finalized name at the end
		self._oldNames = defaultdict(list)	# full-length Allele Code before naming begins (used to track changed codes)
		self._treeHasBeenBuilt = False		# check for tree building status
	
	@staticmethod
	def CDCName(part):
		"""
		CDCName:  prepend appropriate suffix and version to name list converted to string
			Arguments:
				part:	list[int] --> partial Allele Code
			Returns:
				string (concatenated input Allele Code (part) prepended with user-supplied prefix and built-in version)
		"""
		# verify input is a list
		if not isinstance(part, list):
			raise AssertionError('Need list to create code!')
		# add prefix and version onto dot-separated name and return it
		code = '{}{} - {}'.format(prefix, version, Tree.NameToStr(part))
		return code

	def FinalizeName(self, key, name):
		"""
		FinalizeName:  set input key's name (list) to input name (list)
			Arguments:
				key:	string
				name:	list[int]
		"""
		# get node corresponding to input name list
		node = self.Traverse(name)
		# verify input key has been named
		assert(node.IsNamed(key))
		# set the key's name to input name in self._names
		self._names[key] = name

	def FinalizeCDCNames(self):
		"""
		FinalizeCDCNames:	for each NamedNode in tree, set its Keys' name (Allele Code) to the level containing more than one Key
			Returns:
				iterator (tuple of key, Allele Code of appropriate length, and whether or not it's a full-length code)
		"""
		# for each child Node of head Node,
		for node in self._tree.GetChildrenNodes():
			# for each NamedNode of head Node's children,
			for namedNode in node.DFSNamed():
				# get its full-length Allele Code as a list
				cdcname = namedNode.NTraverse()
				part = 0
				# if a digit is -1, then set 'part' to the previous position, so it's truncated below
				for i, val in enumerate(cdcname, 1):
					if val > 0:
						part = len(cdcname) - i
						break
				# yield tuple of each key in NamedNode with [possily] truncated name and whether or not the code is full-length
				for key in namedNode.GetNamedChildren():
					partialName = Tree.CDCName(self.GetPart(key, part))
					complete = part == Tree.DEPTH
					yield (key, partialName, complete)

	def GetName(self, key):
		"""
		GetName:  return Allele Code for input key as list of ints, or empty list if not in self._names dict
			Arguments:
				key:  string
			Returns:
				list (Allele Code for input key as list of ints)
		"""
		return self._names.setdefault(key, [])

	def GetNames(self):
		"""
		GetNames:  return self._names (dict)
			Returns:
				dict ({key:[int, int, int, ...]...})
		"""
		return self._names

	def GetPart(self, key, level):
		"""
		GetPart:  return input key's Allele Code as list up to input digit number (level)
			Arguments:
				key:	string
				level:	int
			Returns:
				list
		"""
		return self.GetName(key)[:level]

	def GetStrName(self, key):
		"""
		GetStrName:  return input key's Allele Code as string, using NameToStr function
			Arguments:
				key:  string
			Returns:
				string
		"""
		return Tree.NameToStr(self.GetName(key))

	def HasName(self, key):
		"""
		HasName:  return whether or not input key is in self._names dict
			Arguments:
				key:  string
			Returns:
				bool
		"""
		return key in self._names

	def Load(self, flobj):
		"""
		Load:  get json-formatted Tree (dict) from file and set initialized variables (__init__) to corresponding values
			Arguments:
				flobj: file object --> open file object for reading
		"""
		# read file into dictionary and verify type
		database = json.load(flobj)
		assert isinstance(database, dict)
		# construct self._tree using Node.Load constructor
		self._tree = Node.Load(None, database['tree'], self.DEPTH)
		# update member variables with tree file contents
		self._treeHasBeenBuilt=True
		self._names.update(database['names'])
		self._oldNames.update(database['names'])
		
	@staticmethod
	def NameToStr(parts):
		"""
		NameToStr:  return dot-separated string comprised of input digit list of ints
			Arguments:
				parts:	list[int] --> list of integers that may or may not be a full-length Allele Code
			Returns:
				string (dot-separated string of concatenated input integer list items)
		"""
		# make sure input is an iterable (list or tuple)
		assert(isinstance(parts, list) or isinstance(parts, tuple))
		# return concatenated digits separated by dots ('.')
		return '.'.join(list(map(str, parts)))

	@staticmethod
	def NameFromStr(name):
		"""
		NameFromStr:  convert into name string into list of ints
			Arguments:
				name:  string
			Returns:
				list (int-converted list split from input string on dots ('.'))
		"""
		# make sure input is indeed a string
		try:
			assert(isinstance(name, str))
		except AssertionError:
			print(name)
		
		# return int-converted values from input split on dots('.')
		return list(map(int, name.split('.')))
		
	@staticmethod
	def CDCNameFromStr(name):
		"""
		CDCNameFromStr:	use NameFromStr function above to return input name string as a list of ints
						after verifying its prefix is appropriate
			Arguments:
				name:  string
			Returns:
				list of ints
		"""
		# make sure input name is indeed a string
		try:
			assert(isinstance(name, str))
		except AssertionError:
			print(name)
		
		# verify input prefix+version prependage is what it should be
		codePrefix = '{}{} - '.format(prefix, version)
		if name[:len(codePrefix)]!=codePrefix:
			# return empty list if not
			return []
		
		# pass input name without prefix to NameFromStr function, and return the list[int] that comes back
		return Tree.NameFromStr(name[len(codePrefix):])

	def UpdateNamed(self, key, newKey):
		"""
		UpdateNamed:  	replace key with newKey in all Nodes, starting with NamedNode at the bottom,
						working up the tree to replace matching "Preferred" keys if applicable
			Arguments:
				key:	 string --> old key to replace
				newKey:  string --> key to replace the old one with
		"""
		# get full-length Allele Code of input key as list
		location = self._names.get(key, None)
		if location is not None:
			if self._treeHasBeenBuilt:
				# get the terminal NamedNode
				named_node = self.Traverse(location)
				if not isinstance(named_node, NamedNode):
					# raise an error if we got to the last Node and it's not a NamedNode class
					raise RuntimeError('Have partial key for which logic requested removal')
				# run duplicate function on terminal NamedNode
				named_node.UpdateNamed(key, newKey)
				
				# then replace the key as "Preferred" on the way back up the tree
				parent = self._tree
				for ID in location:
					# get Node for each digit in code
					child = parent.GetChild(ID)
					# raise error if not found
					if not child:
						raise RuntimeError('Have missing internal node along the path of a named key'
						' for which logic requested removal')
					# otherwise, replace its "Preferred" status if applicable
					if child.Preferred==key:
						child.Preferred = newKey
					# update parent Node variable so next digit's Node is found
					parent = child
					
			# Lastly, replace the key in the names dict
			self._names[newKey] = self._names[key]
			del self._names[key]
	
	def Save(self, flobj):
		"""
		Save:  write {'names':{}, 'tree':{}} dict into input open file object
			Arguments:
				flobj:  file object --> open stream to write data to file
		"""
		json.dump({ 'names': self._names, 'tree': self._tree.Save() }, flobj)

	def Traverse(self, ids):
		"""
		Traverse:  returns head Node's Traverse function results using input Node IDs (ids)
			Arguments:
				ids:  list[string] --> list of Node IDs to pass to head Node
			Returns:
				Node or NamedNode (terminal child Node corresponding to last digit of ids list)
		"""
		assert(isinstance(ids, tuple) or isinstance(ids, list))
		return self._tree.Traverse(ids)

	def Tree(self):
		"""
		Tree:  self-returning member function
			Returns:
				Tree
		"""
		return self._tree

	def __len__(self):
		"""
		__len__:  returns total named keys in Tree if len() function is applied to Tree variable
			Returns:
				int
		"""
		return len(self._names)


class Node(object):
	"""
	Node:  base class for elements managed in Tree class to hold hierarchical naming structure
	"""
	def __init__(self, ID, level, parent):
		"""
		Initialize with input variables and empty/0/None objects
		"""
		self._ID = ID			# numerical ID of node
		self._level = level		# digit/position in period-separated code
		self._parent = parent	# parent node
		self._children = {}		# child nodes
		self._diameter = 0		# maximum weighted distance between profiles in node
		self._preferred = None	# founder/initiator Key
		
		#for caching purposes only!
		self._entryKeys = None
		
	@property
	def Diameter(self):
		"""
		Diameter:  return Node's max distance from self._diameter string
			Returns:
				float
		"""
		return self._diameter	
	
	@Diameter.setter
	def Diameter(self, diameter):
		"""
		Diameter.setter:  set Node diameter (i.e. max within-Node distance) from input float
			Arguments:
				diameter:  float
		"""
		self._diameter = diameter
		
	@property
	def Preferred(self):
		"""
		Preferred:  return current self._preferred key (i.e. Node founder)
			Returns:
				string (self._preferred)
		"""
		return self._preferred	
	
	@Preferred.setter
	def Preferred(self, preferred):
		"""
		Preferred.setter:  set Preferred member variable from input "preferred" string
			Arguments:
				preferred:  string --> founder key for Node
		"""
		self._preferred = preferred
		
	@property
	def EntryKeys(self):
		"""
		EntryKeys:  return self._entryKeys list. if None, compile with iter_entrykeys()
			Returns:
				list (self._entryKeys)
		"""
		if self._entryKeys is None:
			self._entryKeys = list(self.iter_entrykeys())
		return self._entryKeys
		
	@EntryKeys.setter
	def EntryKeys(self, entryKeys):
		"""
		EntryKeys.setter:  set self._entryKeys from input list
			Arguments:
				entryKeys:  list[string] --> list of keys to add to Node's _entryKeys list
		"""
		self._entryKeys = entryKeys
		
	def Save(self):
		"""
		Save: return dictionary of member variables to be written to .json file. calls _Save function in case it's a NamedNode
			Returns:
				dict (dictionary containing current Node variables with child Node variables nested in 'children' dict)
		"""
		# initialize 'info' dict for current Node, recursing for all child Nodes found
		info = {
			'ID': self._ID,
			'level': self._level,
			'diameter': self._diameter,
			'preferred': self._preferred,
			'children': { childId: 
				child.Save() for childId, child in self._children.items()
			}
		}
		# if current Node is a NamedNode, _Save will add 'details' into the info dict, so add it here
		detailedInfo = self._Save()
		if detailedInfo:
			info['details'] = detailedInfo
		
		# return top-down compiled info for Node
		return info
	
	@staticmethod
	def Load(parent, info, depth):
		"""
		Load:  create new Node (GetNode function), and set variables according to input info dict
			Arguments:
				parent:  Node --> parent Node to use for GetNode function
				info:	 dict --> dictionary with Node-specific Keys for filling new Node's values
				dept:	 int --> integer for GetNode function to create a new Node
			Returns:
				Node (new Node created by GetNode function then filled in here)
		"""
		# create new Node using input args
		node = GetNode(info['ID'], info['level'], parent, depth)
		# set new Node's variables
		node._diameter = info['diameter']
		node._preferred = info['preferred']
		if 'details' in info:
			# if 'details' is a Key in info dict, it's a NamedNode, so run _Load function
			node._Load(info['details'])
		node._children = { int(childId): Node.Load(node, childInfo, depth) for childId, childInfo in info['children'].items() }
		
		# return new, completed Node
		return node
	
	def _Save(self):
		"""
		_Save:  return None here; overloaded in NamedNode subclass for detailed info
		"""
		return None
		
	def _Load(self, info):
		"""
		_Load:  do nothing here; overloaded in NamedNode child class
		"""
		return
		
	def AddChild(self, ID, node):
		"""
		AddChild:  add input node into self._children with input ID if not already there
			Arguments:
				ID:  	int
				node:  	Node
		"""
		# verify input ID isn't already there
		assert self._children.get(ID, None) is None
		# "adopt" input node by setting its parent to self
		node.SetParent(self)
		# add input node to self._children under input ID
		self._children[ID] = node

	def DeleteChild(self, ID):
		"""
		DeleteChild:  removes child Node with input ID
			Arugments:
				ID:  int
		"""
		del self._children[ID]

	def DFSNamed(self):
		"""
		DFSNamed:  yield NamedNodes at or below self Node
			Returns:
				NamedNode as iterator
		"""
		if isinstance(self, NamedNode):
			# if already at a NamedNode, return self and stop
			yield self
		else:
			# otherwise, yield all child nodes that are NamedNodes
			toCheck = list(self.GetChildrenNodes())

			while(len(toCheck)):
				ob = toCheck.pop()
				if isinstance(ob, NamedNode):
					# if this child Node is a NamedNode, yield it
					yield ob
				else:
					# otherwise, add its child Nodes to list to keep it going
					toCheck.extend(ob.GetChildrenNodes())
					
	def iter_entrykeys(self):
		"""
		iter_entrykeys:  yield keys from all NamedNodes at or below self Node
			Returns:
				string (key) as iterator
		"""
		for namedNode in self.DFSNamed():
			for key in namedNode.GetNamedChildren():
				yield key

	def GetChild(self, ID):
		"""
		GetChild:  return child node with input ID
			Arguments:
				ID:  int
			Returns:
				Node
		"""
		return self._children.get(ID, None)

	def GetChildren(self):
		"""
		GetChildren:  return self._children
			Returns:
				dict
		"""
		# return dictionary of child nodes
		return self._children

	def GetChildrenKeys(self):
		"""
		GetChildrenKeys:  return child node IDs
			Returns:
				list of ints
		"""
		return self._children.keys()

	def GetChildrenNodes(self):
		"""
		GetChildrenNodes:  return child node values/nodes
			Returns:
				list of Nodes
		"""
		return self._children.values()

	def ID(self):
		"""
		ID:  returns self._ID
			Returns:
				int
		"""
		return self._ID

	def Level(self):
		"""
		Level:  returns self._level
			Returns:
				int
		"""
		return self._level

	def NewChildNode(self):
		"""
		NewChildNode:  add child node with ID one greater than current max child ID, or 1 if no children present
			Returns:
				Node (new node if current level isn't within 1 of tree depth; new NamedNode otherwise if it is)
		"""
		# determine new Node's ID
		if not len(self._children):
			# if no other children, then new node ID is 1
			nextCluster = 1
		else:
			# otherwise, it's the max ID of existing children plus 1
			nextCluster = max(self._children.keys()) + 1
		
		if self._level == Tree.DEPTH - 1:
			# add as NamedNode if current node is penultimate
			self._children[ nextCluster ] = NamedNode(nextCluster,
														self._level+1, self)
		else:
			# otherwise, add as a regular Node
			self._children[ nextCluster ] = Node(nextCluster, 
													self._level+1, self)
		# return newly created node
		return self.GetChild(nextCluster)

	def MergeNodes(self, nodes):
		"""
		MergeNodes:  combine input nodes into the largest one, renaming each lower node consumed in the process to avoid ID conflicts
			Arguments:
				nodes:  list[Nodes] --> list of Nodes to be combined
			Returns:
				Node (largest Node in list to merge, now containing all other nodes in input list)
		"""
		allNodes = [self.GetChild(i) for i in nodes]
		minObjs = {}
		# find largest node by total keys
		for node in allNodes:
			size = sum(n.TotalChildCount() for n in node.DFSNamed())
			minObjs[node] = size
		maxObj = max(minObjs, key=lambda x: minObjs[x])
		
		# remove largest node from working dict
		del minObjs[maxObj]

		if self._level == Tree.DEPTH - 1:
			#all nodes at this level are named nodes, so we need to do a special merge
			for node in minObjs:
				# add all Keys from other nodes into largest one...
				for key in node.GetNamedChildren():
					maxObj.AddNamedChild(key)
				# before deleting those nodes from tree
				self.DeleteChild(node.ID())
		else:
			# rename smaller nodes to next largest integers after largest node's ID
			toStart = max(maxObj.GetChildrenKeys()) + 1
			for node in minObjs:
				for i, child in enumerate(node.GetChildrenNodes(), toStart):
					child.SetID(i)
					maxObj.AddChild(i, child)

				toStart = i + 1
				# remove old node from tree, since everything's been copied to largest node
				self.DeleteChild(node.ID())
		# update Allele Codes for all keys in combined Node
		for node in maxObj.DFSNamed():
			wgst = node.ValidateNames()
			for key in node.GetNamedChildren():
				Tree.CURRENT.FinalizeName(key, wgst)
				
		# return Node resulting from merge
		return maxObj

	def SetID(self, ID):
		"""
		SetID:  update self._ID with input ID if it's an integer
			Arguments:
				ID:  int --> the value to set self._ID to
		"""
		assert isinstance(ID, int)
		self._ID = ID

	def SetParent(self, parent):
		"""
		SetParent:  set Node's parent if input parent is a Node object
			Arguments:
				parent:  Node --> Node to set current Node's self._parent to
		"""
		assert isinstance(parent, Node)
		self._parent = parent

	def TotalChildCount(self):
		"""
		TotalChildCount:  return total child nodes
			Returns:
				int
		"""
		return len(self._children) + self.NamedChildCount()
		
	def NamedChildCount(self):
		"""
		NamedChildCount:  return 0, since this is not a NamedNode and that subclass has a duplicate function
			Returns:
				int
		"""
		return 0

	def Traverse(self, ids):
		"""
		Traverse:  return terminal child node associated with last digit of input code
			Arguments:
				ids: list[ints] --> list of Node IDs
			Returns:
				Node (None if none exist)
		"""
		if len(ids) == 1:
			return self.GetChild(ids[0])
		
		else:
			child = self.GetChild(ids[0])
			
			if child:
				return child.Traverse(ids[1:])
			else:
				return None

	def NTraverse(self):
		"""
		NTraverse:  return node IDs as list from child node with -1 as ID if no child nodes present
			Returns:
				list[ints]
		"""
		if self._parent == None:
			# parent Node doesn't exist, so
			if self.TotalChildCount() == 1:
				# return [-1] if head node of entire tree
				return [-1]
			else:
				# return Node ID otherwise
				return [self.ID()]
		else:
			if self.TotalChildCount() == 1:
				return [-1] + self._parent.NTraverse()
			else:
				return [ self.ID() ] + self._parent.NTraverse()

	def RTraverse(self):
		"""
		RTraverse:  return node IDs as list from child node up until no parent exists
			Returns:
				list[ints]
		"""
		if self._parent == None:
			return []
		else:
			return [self._ID] + self._parent.RTraverse()

	def __repr__(self):
		"""
		__repr__:  return shorthand for Node if printed directly
			Returns:
				string
		"""
		return ('Class: {} | ID: {} | Level: {} | Parent: {}'.format(type(self),
																	self._ID,
																	self._level,
																	self._parent.ID()))


class NamedNode(Node):
	"""
	NamedNode:  class inheriting Node elements with added _namedChildren dict to hold corresponding keys with full-length code, 
				should more than one key be present here
	"""
	
	def __init__(self, ID, level, parent):
		"""
		Initialize class using parent class Node constructor followed by additional variables
		"""
		super(self.__class__, self).__init__(ID, level, parent)	# Node constructor creating ID and level elements with parent input Node
		self._wgst = self.RTraverse()[-1::-1]					# list of ints, representing the full-length Allele Code at this node
		self._namedChildren = {}								# empty dict to hold keys contained in this Node
		
	def _Save(self):
		"""
		_Save: return self._namedChildren as key:list containing all keys in node
			Returns:
				dict --> key = 'namedChildren', values = list of node's keys
		"""
		return { 'namedChildren': list(self._namedChildren.keys()) }
		
	def _Load(self, info):
		"""
		_Load:  set self._namedChildren equal to {key:self._wgst...} for every key in input info dict
			Arguments:
				info:  dict --> json "info" element of saved tree file
		"""
		self._namedChildren = { key: self._wgst for key in info['namedChildren'] }
		
	def AddNamedChild(self, key):
		"""
		AddNamedChild: add input key to node with node's code and return that code
			Arguments:
				key:  string --> key to add
			Returns:
				list of ints
		"""
		self._namedChildren[key] = self._wgst
		return self.GetChildName(key)
	
	def GetChildName(self, key):
		"""
		GetChildName:  return code for input key if present (None otherwise)
			Returns:
				list of ints
		"""
		return self._namedChildren.get(key, None)

	def GetNamedChildren(self):
		"""
		GetNamedChildren:  return self._namedChildren dictionary
			Returns:
				dict
		"""
		return self._namedChildren
		
	def NamedChildCount(self):
		"""
		NamedChildCount:  return length of self._namedChildren
			Returns:
				int
		"""
		return len(self._namedChildren)
	
	def IsNamed(self, key):
		"""
		IsNamed:  return whether input key is in self._namedChildren dictionary
			Arguments:
				key:  string --> dictionary key to search for
			Returns:
				bool
		"""
		return key in self._namedChildren

	def UpdateNamed(self, key, newKey):
		"""
		UpdateNamed:  replace input key with newKey in self._namedChildren
			Arguments:
				key:  	string --> key to replace
				newKey: string --> replacement key
		"""
		# add newKey to self._namedChildren in old key's placed
		self._namedChildren[newKey] = self._namedChildren[key]
		# and delete old key from self._namedChildren
		del self._namedChildren[ key ]

	def ValidateNames(self):
		"""
		ValidateNames:  fill self._wgst list with corresponding digits for each index
			Returns:
				list of ints
		"""
		for i, val in enumerate(self.RTraverse()[-1::-1]):
			self._wgst[i] = val

		return self._wgst

	def Address(self):
		"""
		Address:  return node's code
		"""
		return self._wgst


def GetNode(ID, level, parent, depth):
	"""
	GetNode:  make and return new Node with input values, or NamedNode if at final digit / maximum depth
		Arguments:
			level:  int --> corresponding index of Allele Code (e.g. 1 = 1st digit, 2 = 2nd digit, etc.)
			parent:	Node --> parent node to which this new one is a child
			depth:	int --> maximum number of digits in Allele Code (i.e. len(thresholds[prefix]))
		Returns:
			Node (if level!=depth) or NamedNode (if level==depth)
	"""
	if level==depth:
		return NamedNode(ID, level, parent)
	else:
		return Node(ID, level, parent)


class AlleleCalls(object):
	"""
	AlleleCalls:  class to hold allele profiles with member functions to deal with adding, removing, and saving to file
	"""
	def __init__(self):
		"""
		Initialize class with empty variables
		"""
		self._alleleCalls = {}	# key:profile dict holding allele calls in order specified by coreLoci[] list at top
		self._index = {}		# key:# dict, where # indicates the matrix.#.gzip file holding the corresponding allele profile
		self._matrices = []		# profiles loaded from various matrix.#.gzip files
		
	def _Convert(self):
		"""
		_Convert:	write self._alleleCalls to respective matrix.#.gzip file while more than 1000 profiles remain in
					self._alleleCalls.  Leftovers will be written to calls.gzip by Save function.
		"""
		# total number of profiles to write to file
		BLOCKSIZE = 1000
		
		while len(self._alleleCalls)>BLOCKSIZE:
			keys = list(self._alleleCalls.keys())
			matrix = {}
			# gather 1000 profiles into matrix dict
			for i, key in enumerate(keys[:BLOCKSIZE]):
				# add allele profile to matrix dict
				matrix[key] = self._alleleCalls[key]
				# and its position in self._matrices to self._index
				self._index[key] = (len(self._matrices), i)
				# then delete from self._alleleCalls
				del self._alleleCalls[key]
				
			# append 1000-key matrix to self._matrices
			self._matrices.append(matrix)
			
			# write 1000-key matrix to matrix.#.gzip file, where # is the current length of self._matrices minus 1
			with gzip.open(os.path.join(self._path, 'matrix.{}.gzip'.format(len(self._matrices)-1)), 'wb') as flobj:
				flobj.write(str(matrix).replace("'", '"').encode())
				
	def Load(self, path):
		"""
		Load:  reads calls.gzip and index.gzip from 'current' folder into self._alleleCalls and self._index, respectively
			Arguments:
				path:  string --> path to 'current' folder in data directory
		"""
		# update self._path to input path
		self._path = path
		
		# read calls.gzip into content dict,
		flName = os.path.join(path, 'calls.gzip')
		if os.path.exists(flName):
			with gzip.open(flName, 'rb') as flobj:
				content = json.load(flobj)
				# then convert allele profiles to list of ints when loading into self._alleleCalls
				for key, profile in content.items():
					self._alleleCalls[key] = list(map(int, profile))
		
		# read index.gzip into self._index
		flName = os.path.join(path, 'index.gzip')
		if os.path.exists(flName):
			with gzip.open(flName, 'rb') as flobj:
				self._index = json.load(flobj)
				
		self._matrices = []
	
	def _LoadMatrix(self, i):
		"""
		_LoadMatrix: read matrix.[i].gzip into self._matrices[i]
			Arguments:
				i:  integer --> name in the middle of saved matrix.[i].gzip file, holding saved allele calls
		"""
		# make self._matrices at least as long as input i
		while len(self._matrices) <=i:
			self._matrices.append(None)
			
		# read matrix.i.gzip ito self._matrices[i]
		with gzip.open(os.path.join(self._path, 'matrix.{}.gzip'.format(i)), 'rb') as flobj:
			self._matrices[i] = json.load(flobj)
	
	def Save(self, path):
		"""
		Save:  save Allele Calls to file with keys mapped to file names in index file
			Arguments:
				path:  string --> path to 'current' folder for allele calls
		"""
		# update self._path to input path
		self._path  = path
		
		# save bulk of calls to matrix.#.gzip files in data directory
		# (matrix.#.gzip files created after 1000 keys have been named)
		self._Convert()
		
		# convert leftovers to lists
		for key, profile in self._alleleCalls.items():
			self._alleleCalls[key] = list(profile)
		
		# and write to calls.gzip file, converting self._alleleCalls dict to byte-encoded, jason-serializable string
		with gzip.open(os.path.join(path, 'calls.gzip'), 'wb') as flobj:
			flobj.write(str(self._alleleCalls).replace("'", '"').encode())
		
		# save the index to file too, converting as with self._alleleCalls above
		with gzip.open(os.path.join(path, 'index.gzip'), 'wb') as flobj:
			flobj.write(str(self._index).replace("'", '"').encode())
			
	def Add(self, key, calls):
		"""
		Add:  adds input allele calls (calls) to self._alleleCalls if not already in self._alleleCalls or self._index
			Arguments:
				key:  string --> key for self._alleleCalls and self._index dictionaries
				calls:	list --> allele profile to add to self._alleleCalls
		"""
		# Log if input key is already part of the class
		assert(key not in self._alleleCalls), "Hmm, are you trying to update a key here?"
		assert(key not in self._index), "Hmm, are you trying to update a key here?"
		# add input allele profile to self._alleleCalls for input key
		self._alleleCalls[ key ] = calls

	def GetCalls(self, key):
		"""
		GetCalls:  return allele calls for input key from matrix, if present; None otherwise
			Arguments:
				key:  string
			Returns:
				list (allele calls) if already loaded or loads it first then returns it if not
		"""
		idx = self._index.get(key, None)
		if idx:
			if idx[0] >= len(self._matrices) or self._matrices[idx[0]] is None:
				# if input key's allele profile isn't loaded yet, load the matrix file corresponding to it
				self._LoadMatrix(idx[0])
				
			return self._matrices[idx[0]][idx[1]]
		
		return self._alleleCalls.get(key, None)

	def HasKey(self, key):
		"""
		HasKey:  return whether input key is in either _index or _alleleCalls dictionaries
			Arguments:
				key:  string --> key to look for in _index and _alleleCalls dictionaries
			Returns:
				bool
		"""
		return key in self._index or key in self._alleleCalls
		
	def __len__(self):
		return len(self._alleleCalls) + len(self._index)


def GetDistance(p1, p2):
	"""
	GetDistance: return total differences between profile 1 (p1) and profile 2 (p2)
					weighted by total indexes called in both (> 0)
		Arguments:
			p1:  list[int] --> allele profile 1
			p2:  list[int] --> allele profile 2
			
		Returns:
			float
	"""
	# update total distances calculated
	global cntDistancesCalculated
	cntDistancesCalculated += 1
	# compute total loci called in both profiles
	nCommon = sum([p1[i]>0 and p2[i]>0 for i in range(len(p1))])
	# calculate total allele calls that differ between profiles
	nDiff = sum([p1[i] != p2[i] for i in range(len(p1)) if p1[i] and p2[i]])
	# if profiles share called loci (i.e. nCommon > 0),
	if nCommon:
		# return percent different
		return 100.0 * float(nDiff) / float(nCommon)
	else:
		# if not, return 100% different
		return 100.


#=========================== NAMING FUNCTION ================================#
def CalcName(named, tree, alleles, unNamedEntry, thresholds, corePercent):
	"""
	CalcName: Primary naming function.  Determines if Key to be named (unNamedEntry) belongs in any existing
				Node of named keys (named) in input Tree (tree), using input distance thresholds (thresholds),
				allele profiles (alleles), and minimum acceptable percent called loci (corePercent).
		Arguments:
			named:  		dict --> {key:[code]...} object of everything already named in input tree
			tree:			Tree --> tree class to use for assigning Allele Codes
			alleles:		dict --> {key:[profiles]...} object of both named and unnamed allele profiles
			unNamedEntry:	string --> key of new allele profile to assign an Allele Code
			thresholds:		list --> list of numbers (float) for iterative distance calculation
			corePercent:	float --> total percent of input allele profile that must be > 0
			
		Returns:
			updated tree (Tree) with new allele profile added at appropriate Node
	"""
	# Create a local copy of the named keys
	namedEntries = named if isinstance(named, list) else list(named)

	# keep distance locally, only calculate those you really need
	distances = {}
	index = { namedEntryKey: i for i, namedEntryKey in enumerate(namedEntries) }
	
	def GetDistanceToNamedEntry(unNamedEntryKey, namedEntryKey):
		"""
		GetDistanceToNamedEntry:	helper function to add total allelic profile differences to local distances dict
									for input unnamed and named Keys
			Arguments:
				unNamedEntryKey:	string --> Key of allele profile that HAS an Allele Code
				namedEntryKey:		string --> Key of allele profile undergoing Allele Code assignment
				
			Returns:
				float --> result of GetDistance function between allele profiles corresponding to input unnamed and named Keys
		"""
		# add unNamedEntryKey to distances dict if not already there with -1 for all pairwise distances
		if unNamedEntryKey not in distances:
			distances[unNamedEntryKey] = [-1.0]*len(namedEntries)
		
		# myDistances: shorthand for input unNamedEntryKey's distances (all -1's at start)
		myDistances = distances[unNamedEntryKey]
		
		# if distance between input Keys is already determined, return it
		i = index[namedEntryKey]
		if myDistances[i]>=0: 
			return myDistances[i]
		
		# otherwise, update local distance dict with GetDistance result and return it
		v1 = alleles.GetCalls(namedEntryKey)
		v2 = alleles.GetCalls(unNamedEntryKey)
		d = GetDistance(v1, v2)
		
		myDistances[i] = d
		
		return d

	def IsInCluster(unNamedEntryKey, node, threshold, corePercent):
		"""
		IsInCluster:	returns whether input unNamedEntryKey belongs in current node
						
			Arguments:
				unNamedEntryKey:  string --> Key of allelic profile being assessed
				node:	Node object --> node containing keys to be added or excluded by comparing distance to input threshold
				threshold:  float --> maximum distance value allowed for adding to input Node
				corePercent:  float --> percent of allelic profile that must be > 0 to pass QC
			Returns:
				bool -->	1. True if distance to founder ("Preferred") is less than input threshold
							2. False if 2*node diameter plus 10% buffer is greater than input threshold
							3. True if in the "grey zone", where unNamedEntryKey's profile is more distant
								from founder than input threshold but less than threshold for another member of node
							4. False if all of the above fails
		"""
		# distance to preferred cluster sample
		d = GetDistanceToNamedEntry(unNamedEntryKey, node.Preferred)
		if d <= threshold: 
			return True
		
		# if the sample is too far off, just say it's too far off
		if d - node.Diameter - 2.0*(100.0 - 100.0*corePercent) > threshold:
			# calculated distance is more than Node diameter plus 10% buffer
			return False
		else:
			# calculated distance might be in cloud around preferred Key, so compare to everything in the Node
			for sample2 in node.EntryKeys:
				if GetDistanceToNamedEntry(unNamedEntryKey, sample2) <= threshold:
					# if distance to any existing Keys is under threshold, then it belongs there
					return True
			# if we haven't found a single match, then it doesn't belong
			return False
			
	def GetMaxDistanceToNode(entryKey, node):
		"""
		GetMaxDistanceToNode:  return maximum distance of input key's (entryKey) allele profile to all others in the input node
			Arguments:
				entryKey:	string --> Key of unnamed profile
				node:		Node --> Node holding all allelic profiles to compare to
			Returns:
				float --> maximum distance of one-vs-all comparisons between input key's allele profile and all others in input node
		"""
		return max(GetDistanceToNamedEntry(entryKey, m) for m in node.EntryKeys)
		
	# Make sure thresholds are sorted biggest first
	thresholds.sort(key=lambda x: -x)

	# headNode = top node of tree (all functional nodes are children of headNode)
	headNode = tree.Tree()

	# Ensure unNamedEntryKey has no name
	patternName = tree.GetName(unNamedEntry)
	assert len(patternName) == 0

	# currentNode = placeholder for headNode as unNamedEntryKey is placed into child nodes
	currentNode = headNode

	# See if unnamed entry belongs to any existing Nodes at each threshold
	clusterIdAtPreviousLevel = ''
	for level, threshold in enumerate(thresholds, 1):

		# Create the clusters at this level:
		closestClusters = set()
		matchingNodes = []

		for node in currentNode.GetChildrenNodes():
			# if unNamedEntry key really belongs in this node,
			if IsInCluster(unNamedEntry, node, threshold, corePercent):
				# add digit to closestClusters set
				keyToGet = node.Preferred
				partial = tuple(tree.GetPart(keyToGet, level))								
				closestClusters.add(partial)
				# and node itself to list of verified matches
				matchingNodes.append(node)
						
		# if only one matching cluster/node,
		if len(closestClusters) == 1:
			# add node's Allele Code digit to patternName list,
			pattern = closestClusters.pop()
			patternName.append(pattern[-1])
			# update working node (currentNode) to this matched node, so later iteration continue down this branch
			# (hierarchical naming)
			currentNode = matchingNodes[0]
			# update the node diameter to the max distance between its current value 
			# and all pairwise comparisons with the new unNamedEntry,
			currentNode.Diameter = max(currentNode.Diameter, GetDistanceToNamedEntry(unNamedEntry, currentNode.Preferred))
			# and add unNamedEntry key to updated Node's list of keys
			currentNode.EntryKeys.append(unNamedEntry)
						
		# if NO matching nodes were found,
		elif len(closestClusters) == 0:
			# then make a new one
			currentNode = currentNode.NewChildNode()
			# add new node's ID to patternName list for unNamedEntry
			patternName.append(currentNode.ID())
			# set initial node diameter to 0, since it contains only a single allele profile
			currentNode.Diameter = 0
			# set the found ("Preferred") to this unNamedEntry key
			currentNode.Preferred = unNamedEntry
			# and add that key to the node's list of keys
			currentNode.EntryKeys = [unNamedEntry]		

		# if unNamedEntry matched more than one node, then all of them need to be merged
		# into a new one (single-linkage principle)
		else:
			# find largest matching node by total keys within it
			maxSize = max(len(node.EntryKeys) for node in matchingNodes)
			# get index of largest node in that list
			nodePickIdx = next(i for i, node in enumerate(matchingNodes) if len(node.EntryKeys)==maxSize)
			# get largest node's "Preferred" key
			newPreferred = matchingNodes[nodePickIdx].Preferred
			# and finally get its diamter
			diameter = matchingNodes[nodePickIdx].Diameter			
			
			# update it with new unNamedEntry
			diameter  = max(diameter, GetDistanceToNamedEntry(unNamedEntry, newPreferred))
			
			# and then all the other nodes			
			for i, node in enumerate(matchingNodes):
				# skipt largest node, since that's what we're working around
				if i==nodePickIdx: continue
				diameter = max(diameter, GetMaxDistanceToNode(newPreferred, node))				
				
			# collect keys from all matching nodes
			entryKeys = [unNamedEntry]
			for node in matchingNodes:				
				entryKeys.extend(node.EntryKeys)

			# compile terminal digits for matching nodes' Allele Codes
			toMerge =[ c[-1] for c in closestClusters ]
			# reset working node (currentNode) to the results of MergeNodes with terminal digits as input
			currentNode = currentNode.MergeNodes(toMerge)
			# and add the new node's ID to unNamedEntry's Allele Code digit list
			patternName.append(currentNode.ID())
			
			# Update new, merged node using these values
			currentNode.Preferred = newPreferred
			currentNode.Diameter = diameter
			currentNode.EntryKeys = entryKeys
		
	# alert if the above loop completed assessing all thresholds but currentNode didn't end up being a terminal, NamedNode
	assert isinstance(currentNode, NamedNode), "Hmm, we are at the end here but we don't have a NamedNode?"
	
	# add unNamedEntry to "namedChildren" list
	currentNode.AddNamedChild(unNamedEntry)	

	# alert if unNamedEntry's Allele Code wasn't updated to the compiled patternName list
	assert currentNode.GetChildName(unNamedEntry) == patternName,  "Hmm, {} seems to be in the wrong named node? {}  != {} ".format(unNamedEntry, currentNode.GetChildName(unNamedEntry), patternName)

	# in the end, return the edited input tree
	return tree


class Calculator(object):
	"""
	Calculator class:  loads data, performs QC, and functions as a wrapper to call the primary naming function (CalcName) on each new allele profile
	"""
	def __init__(self, runtimeArgs):
		"""
		Create Calculator class with various member variables derived from input runtimeArgs and global variables
		"""
		self._minPres = minpres									# global minpres variable
		self._thresholds = list(thresholds[prefix])				# prefix-specific list of distance thresholds
		self._datadir = runtimeArgs['dirpath']					# path to data directory
		self._logdir = runtimeArgs['logpath']					# path to log diretory
		self._treepath = runtimeArgs['tree']					# path to tree file in 'current' folder
		self._callspath = runtimeArgs['allele calls']			# path to 'current' allele calls folder
		self._changeLogPath = runtimeArgs['changeLog']			# path to today's change_log file
		self._XcodesPath = runtimeArgs['Xcodes']				# path to Xcodes.tsv file, if available
		self._nosave = nosave									# bool:  whether or not to overwrite saved data files
		self._tree = Tree(len(self._thresholds))				# initial Tree class as empty Node
		self._alleleCalls = AlleleCalls()						# empty AlleleCalls class
		self._newAlleleCallsFile = runtimeArgs['new alleles']	# path to file containing new allele profiles for naming
		self._newKeys = set()									# set to track new profiles being named
		self._changedKeys = set()								# set to track previously named profiles whose names change this run
		self._belowQC = defaultdict(list)						# dict to track new profiles that don't pass QC
		self._numChanged = 0									# total codes that have changed this run
		self._output = runtimeArgs['outputFilePath']			# path to output file, where results will be written
		
	def DoRefresh(self, args):
		"""
		DoRefresh:  delete lock file, since run is done, and let user know it's done
		"""
		numChanged = self._numChanged
		
		try:
			lockfile_path = os.path.join(DATA_DIR,
										'{}_nomenclature_srcfiles'.format(prefix),
										'nomenclature.lock')
			os.remove(lockfile_path)
		except:
			pass
			
		# Let user know assignment is done
		if self._nosave:
			print('Assignment done without saving. ({} total code changes)'.format(numChanged))
		else:
			print('Assignment complete with {} changed allele codes.'.format(numChanged))

	def RunCalc(self, args):
		"""
		RunCalc: calls DoCalc function to assign Allele Codes to new profiles, and removes lock file if run without errors
		"""
		global Logger
		Logger = logging.getLogger(__name__)
		
		# Log that run is beginning
		log_message("=== START OF NOMENCLATURE RUN ===")
		
		try:
			lockfile_path = os.path.join(DATA_DIR,
										'{}_nomenclature_srcfiles'.format(prefix),
										'nomenclature.lock')
			# Run DoCalc function below, which does all the work
			self.DoCalc(args)

		except:
			# If error occured in DoCalc, log that it happened and alert user
			log_exception('There was an error! Insertion point=RunCalc')
			if not verbose:
				print('There was an error! Insertion point=RunCalc')

		else:
			# If everything went without error, remove the lock file
			os.remove(lockfile_path)

		finally:
			# No matter what, log that the run has ended with total distances calculated during its execution
			log_message("=== END OF NOMENCLATURE RUN ===")
			global cntDistancesCalculated
			log_message(" --> number of distances calculated in this run: {}".format(cntDistancesCalculated))
			close_handlers(Logger.handlers)

	def WriteResults(self, outDelim):
		"""
		Write results to output file path (self._output), using input delimiter to determine file extension
		"""
		with open(self._output, 'w') as f:
			# start with header line
			f.write('{}\n'.format(outDelim.join(['Key', 'Allele_code'])))
			
			for key, value, complete in self._tree.FinalizeCDCNames():
				# new and changed codes
				if key in self._newKeys or key in self._changedKeys:
					f.write('{}\n'.format(outDelim.join([key, value])))
				# new profiles that didn't pass QC
				if key in self._belowQC:
					criteria = ', '.join(self._belowQC[key])
					f.write('{}\n'.format(outDelim.join([key, 'FAILED QC: {}'.format(criteria) ])))
			
	def DoCalc(self, args):
		"""
		DoCalc: performs Allele Code assignment from user-provided arguments loaded during class creation
		
		Steps:
			1.	Load nomenclature tree and allele profiles from files
			2.	Load new allele profiles from user-provided path
			3.  Perform QC on each profile, and if it passes,
					3a.	Calculate distance to each founder ("preferred" key) profile in nomenclature tree.
					3b.	Create new node and stop if distance is too great, or recursively add to node if below organism-specific threshold(s)
					3c.	Repeat 3a-3b until all thresholds have been assessed
			4.	Save results if --nosave flag was not provided
		"""
		global numChanged
		
		def CheckCore(calls):
			"""
			CheckCore: return None if more 0's than allowed in allele profile or the input allele profile if not
			"""
			if round(float(len(calls) - calls.count(0)) / float(len(calls)),2) < self._minPres:
				return None
			else:
				return calls

		def SaveTree():
			"""
			SaveTree: overwrite .json file at self._treepath with most recent nomenclature results
				Returns:  string --> path to newly written nomenclature tree file (.json format)
			"""
			FILE_PATH = None
			if not os.path.isdir(self._treepath):
				# if tree path isn't just a folder, make new time-stamped file in that directory and delete the old one
				DIR_NAME = os.path.dirname(self._treepath)
				FILE_NAME = 'tree_{}.json'.format(datetime.now().strftime("%Y-%m-%d@%H-%M-%S"))
				FILE_PATH = os.path.join(DIR_NAME, FILE_NAME)

				# use class function to write nomenclature tree to open file
				with open(FILE_PATH, 'w') as f:
					self._tree.Save(f)
				
				# Log successful save and remove old tree from 'current' directory if not an exact match
				# (exact matches occur when less than 1 second has passed between runs)
				log_message('SaveTree: new file created @ {}'.format(FILE_PATH), depth=1)
				if FILE_PATH != self._treepath:
					log_message('SaveTree: removing old tree @ {}'.format(self._treepath), depth=1)
					os.remove(self._treepath)
				
			else:
				# if it is a directory, just write a new file
				FILE_NAME = 'tree_{}.json'.format(datetime.now().strftime("%Y-%m-%d@%H-%M-%S"))
				FILE_PATH = os.path.join(self._treepath, FILE_NAME)

				# Log save operation
				log_message('SaveTree: creating new file @ {}'.format(FILE_PATH), depth=1)
				with open(FILE_PATH, 'w') as f:
					self._tree.Save(f)
					
			# output new tree's file path to reset _treepath member variable after intermittent saves
			return FILE_PATH
					
		
		# Load tree from file if it exists
		if not os.path.isdir(self._treepath):
			log_message('Loading current tree from {}'.format(self._treepath))
			with open(self._treepath, 'rb') as f:
				self._tree.Load(f)
			log_message('Tree loaded successfully from {}'.format(self._treepath), depth=1)
		else:
			# if only a directory, then create a new tree file (should only be for initial naming)
			log_message('No tree file. Creating with SaveTree function.')
			self._treepath = SaveTree()
			# now load as usual
			with open(self._treepath, 'rb') as f:
				self._tree.Load(f)
			log_message('Tree created successfully a {}'.format(self._treepath), depth=1)
			
		# Load allele calls from file storage
		if os.path.isdir(self._callspath):
			log_message('Loading historical allele calls')
			
			self._alleleCalls.Load(self._callspath)
			log_message('Successfully loaded allele calls', depth=1)

		# Set xCodeList from file if applicable
		if self._XcodesPath:
			SetXcodeList(self._XcodesPath)
		if len(xCodeList):
			log_message('{} Xcodes loaded'.format(len(xCodeList)), depth=1)
		else:
			log_message('No Xcodes loaded from file', depth=1)
			
		# Make sure the tree and allele calls file keys match up
		log_message('Checking data integrity')
		if all([alleleKey in self._tree._names.keys() for alleleKey in self._alleleCalls._alleleCalls.keys()]) and \
			all([treeKey in self._alleleCalls._alleleCalls.keys() for treeKey in self._tree._names.keys()]):
			pass
		else:
			# Compile and log keys in one file but not the other
			keysInTreeNotCalls = [alleleKey for alleleKey in self._alleleCalls._alleleCalls.keys() if alleleKey not in self._tree._names.keys()]
			keysInCallsNotTree = [treeKey for treeKey in self._tree._names.keys() if treeKey not in self._alleleCalls._alleleCalls.keys()]
			if len(keysInTreeNotCalls):
				log_error('Keys in tree file missing from calls:', depth=1)
				for key in keysInTreeNotCalls:
					log_error(key, depth=2)
			if len(keysInCallsNotTree):
				log_error('Keys in calls file(s) missing from tree file:', depth=1)
				for key in keysInCallsNotTree:
					log_error(key, depth=2)
					
			raise AssertionError('Error:  Either Keys present in allele file(s) are missing from Tree file or vice versa')
		
		log_message('Successfully checked data integrity', depth=1)
		
		# Run FinalizeCDCNames on current tree, so oldNames attribute will be accurate
		for key, name, complete in self._tree.FinalizeCDCNames():
			self._tree._oldNames[key] = name
		
		# Set algorithm start time
		date = Now()
		
		# Get all the current names
		namedEntries = list(self._tree.GetNames())
		
		# Log start time
		log_message('Beginning calculation at {}'.format(datetime.now().strftime('%I:%M:%S %p')))
		
		# Load new allele profiles in same order as saved data (using coreLoci global variable)
		print('Loading new profiles...')
		newProfiles = LoadProfilesFromFile(newAllelesPath)
		# remove already named Keys
		newProfiles = {key:calls for key, calls in newProfiles.items() if key not in namedEntries}
		numProfiles = len(newProfiles)
		
		# Time to assign names:
		print('Calculating and assigning Allele Codes')
		# progress bar increment, so only 20 dots are printed to terminal
		progressChunk = int(numProfiles/20) if numProfiles > 20 else 1
		nChunks = 0
		
		for i, (key, calls) in enumerate(newProfiles.items()):
			# update progress "bar"
			sys.stdout.write('.'*nChunks + ' '*(20-nChunks) + '{}%'.format(int(float(100*(i+1))/float(numProfiles))) + '\r')
			sys.stdout.flush()
			# update "chunk" variable at every 5% of profiles assessed to keep progress dots at 20 or less
			if i and i % progressChunk == 0:
				nChunks += 1
			
			log_message('Trying to add entry: {}'.format(key), depth=1)

			# Check to see if we have a name in tree file
			if self._tree.HasName(key):
				# Then see if it's missing from allele calls file
				if not self._alleleCalls.HasKey(key):
					# Still in tree file, but not allele calls file,
					# so check QC metrics before adding back to allele calls dictionary
					
					if CheckCore(calls):
						self._alleleCalls.Add(key, calls)
					else:
						log_message('Entry: {} has Allele Code, but is below'
									'the {} presence cutoff'.format(key, self._minPres))

				# Key found in both tree and allele calls files, so skip it
				log_message('This entry has a name, skipping...', depth=1)
				continue
			
			# Perform QC of new allele profile
			log_message('Performing QC...', depth=2)
			
			# add to selection tracker
			self._newKeys.add(key)
			
			qcFail = False
			
			eCalls = CheckCore(calls)

			# Continue if we don't meet QC standards
			if eCalls is None:
				self._belowQC[key].append('CORE')
				qcFail = True
			
			if qcFail:
				log_message('Entry failed QC', depth=3)
				continue
			else:
				self._alleleCalls.Add(key, eCalls)
				log_message('Passed QC', depth=3)
			
			log_message('Assigning name...', depth=2)
			
			# Reset tree with new entry added
			self._tree = CalcName(namedEntries, 
									self._tree, 
									self._alleleCalls, 
									key, 
									self._thresholds, 
									self._minPres)

			# Log successful naming of new profile
			if self._tree.HasName(key):
				namedEntries.append(key)
				log_message('Successfully assigned name!', depth=3)

		# finish by updating progress bar to 100%
		sys.stdout.write('{}{}%\n'.format('.'*22, 100))
		sys.stdout.flush()
		
		# After completing assignment,
		# - save tree and allele calls files if --nosave not provided
		# - save results and changes to output file if -o/--output provided, or print to screen otherwise
		
		# To better tabularize changed codes output
		maxKeyLen = max([0] + [len(key) for key in namedEntries])
		maxOldLen = max([0] + [len(code) for code in self._tree._oldNames.values() if len(code)])
		
		if self._nosave:
			if self._output:
				outDelim = ',' if self._output.endswith('csv') else '\t'
				self.WriteResults(outDelim)
				
			else:
				# output file not given as input, so just print to screen
				for key, value, complete in self._tree.FinalizeCDCNames():
					if key in self._newKeys:
						print('{} {} {}'.format(key, 
												' '*(maxKeyLen - len(key)), 
												value))
					if key in self._changedKeys:
						print('{} {} {} {}\t-->\t{}'.format(key, 
													' '*(maxKeyLen - len(key)), 
													self._alleleCalls._oldNames[key], 
													' '*(maxOldLen - len(currName)), 
													value))
				# follow up with FAILED QC profiles
				for key in self._belowQC:
					if key in self._newKeys:
						print('{} {} {}'.format(key, 
												' '*(maxKeyLen - len(key)), 
												'FAILED QC: {}'.format(', '.join(self._belowQC[key]))))
			
		else:
			# --nosave not provided, so output results to file or terminal
			date = Now()
			
			# Number of names given
			namesGiven = 0
			
			# Number of names given to previously unnamed entries OR those changed by algorithm
			numToSave = 0
			oldNumToSave = 0
			
			# For each key
			for key, name, complete in self._tree.FinalizeCDCNames():
				# Check if the current name is different from previous assignment
				currName = self._tree._oldNames[key]
				if len(currName)==0:
					# new assignments return an empty list, so change currName to string
					# with prefix and version only
					currName = '{}{} - '.format(prefix, version)
				
				#check if the name after applying the blacklisting is different
				xName = CheckXcodeList(name)
				newName = xName
				
				# Log change if one occured
				if '.'.join(currName.split('.')[:len(newName.split('.'))]) != newName:
					# increment numToSave, so incremental save includes changed codes
					numToSave += 1
					
					# log change if not a new assignment
					if len(currName) and not currName.strip().endswith('-'):
						# print changed codes regardless of --verbose flag
						print('{} {} {} {}\t{}\t{}'.format(key, ' '*(maxKeyLen - len(key)), currName, ' '*(maxOldLen - len(currName)), '-->', newName))
						# Add change to date-stamped change_log file with type of change that occured
						changedKeys[key] = {'OldValue':currName, 'NewValue': newName}
						# also indicate type of change
						if currName.endswith('x') or newName.endswith('x'):
							# was previously or is now an Xcode
							changedKeys[key].update({'ChangeType':'X'})
						elif len(currName.split('.')) < len(newName.split('.')):
							# new code is longer,
							if currName == '.'.join(newName.split('.')[:len(currName.split('.'))]):
								# new code is same up to length of old code ("Extended")
								changedKeys[key].update({'ChangeType':'Extended'})
							else:
								# new code changed at some digit ("Merged")
								changedKeys[key].update({'ChangeType':'Merged@{}'.format(min([i for i, c in enumerate(currName.split('.')) if c!=newName.split('.')[i]]))})
						else:
							# other cases:
							cName = currName.split('.')
							nName = newName.split('.')
							if any([cName[i] != nName[i] for i in range(min([len(cName), len(nName)]))]):
								# old code is same length but digit changed ("Merged")
								changedKeys[key].update({'ChangeType':'Merged@{}'.format(min([i for i in range(min([len(cName), len(nName)])) if cName[i]!=nName[i]]))})
							else:
								# another situation, likely removal of a Key from reference files, so code was shortened
								changedKeys[key].update({'ChangeType':'Other'})
					

				if complete:
					namesGiven += 1
					
				# increment numToSave variable, so incremental save includes previously unnamed selected entries
				if key in self._newKeys:
					numToSave += 1
					if not self._output:
						print('{} {} {}'.format(key, 
												' '*(maxKeyLen - len(key)), 
												newName))

				# Save everything after any of the current selection is named and either
				# 1000 of them have been named or previously named entries have changed 
				# to avoid stressing the system and to make sure files are up to date 
				# as best as possible if something breaks
				if numToSave > oldNumToSave and numToSave % 1000 == 0:
					log_message('Intermittent save (numChanged = {})'.format(len(changedKeys)), depth=1)
					print('=====Intermittent save (total changed:  {})====='.format(len(changedKeys)))
					# save over allele profile files
					self._alleleCalls.Save(self._callspath)
					# update tree path to new file to avoid save errors later
					self._treepath = SaveTree()
					# update oldNumToSave, so script won't keep saving files when no selected/changed entries
					# have been encountered while writing to file
					oldNumToSave = numToSave
					
			# print FAILED QC profiles if -o/--output path not provided
			for key in self._belowQC:
				criteria = ', '.join(self._belowQC[key])
				if not self._output:
					print([key, 'FAILED QC: {}'.format(criteria)])

			# One final save after all is done
			self._alleleCalls.Save(self._callspath)
			self._treepath = SaveTree()

			# Write results to output file if given
			if self._output:
				outDelim = ',' if self._output.endswith('csv') else '\t'
				self.WriteResults(outDelim)
			
		# Write changed Allele Codes to file once everything is done
		self._numChanged = len(changedKeys)
		if self._numChanged and len(self._changeLogPath):
			with open(self._changeLogPath, 'a+') as f:
				f.write('{}\n'.format('\n'.join('\t'.join([key, 
															vals['OldValue'], 
															vals['NewValue'], 
															vals['ChangeType']]) for key, vals in changedKeys.items())))
				f.write('=====Assignment Complete ({})====='.format(Now()) + '\n')
				f.close()
				
					
def Run(runtimeArgs):
	"""
	Run:  uses input runtimeArgs to perform Allele Code assignment, using saved tree (.json) and allele calls in data diretory
		Arguments:
			runtimeArgs:  dict --> dictionary with user-defined input parameters and validated directory structure paths
	"""
	# Create Calculator class with input args
	calc = Calculator(runtimeArgs)

	# Do the Allele Code assignment
	calc.RunCalc({})
	
	# Follow up with Refresh action to remove lock file before next run
	calc.DoRefresh({})
			

def Setup():
	"""
	Setup:
		Steps:
			1.	Verifies data directories exist and creates them if not.
			2.	Backs up data files (tree, calls, etc.), then sets RUNTIME_ARGS components to those paths
		
		Returns:
			dict (RUNTIME_ARGS) containing paths to relevant directories (data, logs) and files (tree, calls)
	"""
	
	RUNTIME_ARGS = {
		'dirpath': None,
		'logpath': None,
		'tree': None,
		'allele calls': None,
		'new alleles': newAllelesPath,
		'changeLog': None,
		'nosave': nosave,
		'Xcodes': None,
		'outputFilePath': outputPath if len(outputPath) else None
	}
	
	def ValidateDirectories():
		"""
		ValidateDirectories:	Creates *_nomenclature_srcfiles and *_nomenclature_logs directory structure
								if they don't already exist
		"""
		# set up / verify data directory
		DIR_PATH = os.path.join(DATA_DIR, '{}_nomenclature_srcfiles'.format(prefix))

		if not os.path.isdir(DIR_PATH):
			# if DIR_PATH isn't actually a directory, make one
			log_message('Making srcfiles directory...')
			os.mkdir(DIR_PATH)
			# then fill it with required folders
			for directory in DATA_DIRS:
				path = os.path.join(DIR_PATH, directory)
				os.mkdir(path)

		else:
			# srcfiles directory found, so make subdirectories if missing
			for directory in DATA_DIRS:
				path = os.path.join(DIR_PATH, directory)
				if not os.path.exists(path):
					# Log the missing directory and create it
					log_message('Path: {} did not exist'.format(path))
					os.mkdir(path)

		# update RUNTIME_ARGS with completed data directory path
		RUNTIME_ARGS['dirpath'] = DIR_PATH
		
		# set up / verify log directory
		LOG_PATH = os.path.join(DATA_DIR, '{}_nomenclature_logs'.format(prefix))
			
		if not os.path.isdir(LOG_PATH):
			# as above, if path isn't actually a folder, make it
			log_message('Making logs directory...')
			os.mkdir(LOG_PATH)
			# then fill with required subdirectories
			for directory in LOG_DIRS:
				path = os.path.join(LOG_PATH, directory)
				os.mkdir(path)
		
		else:
			# logs directory found, so make subdirectories if missing
			for directory in LOG_DIRS:
				path = os.path.join(LOG_PATH, directory)
				if not os.path.exists(path):
					# Log the missing directory and create it
					log_message('Path: {} did not exist'.format(path))
					os.mkdir(path)
		
		# update RUNTIME_ARGS with completed log directory path
		RUNTIME_ARGS['logpath'] = LOG_PATH
		
		# Log successful folder creation/validation
		log_message('Successfully validated directories')
		
	def ValidateSrcFiles(RUNTIME_ARGS):
		"""
		ValidateSrcFiles:  Adds file paths to RUNTIME_ARGS dictionary for source files
			Arguments:
				RUNTIME_ARGS:  dict --> original RUNTIME_ARGS dict declared above for editing and later return
			Returns:
				dict (edited RUNTIME_ARGS dict)
		"""
		def Backup(path, hasMultipleFiles):
			"""
			Copies tree file from 'current' subdirectory up one, 
			then returns file path in 'current' subdirectory for use in the algorithm
			"""
			
			# Get the paths ('current' subdirectory)
			PATH_DIR = os.path.join(path, 'current')
			PATH_FILES = [ os.path.join(PATH_DIR, file) for file in \
				os.listdir(PATH_DIR) ]

			# Make sure we only have real files
			PATH_FILES = list(filter(os.path.isfile, PATH_FILES))
			
			# There should never be more than one file in the 'current' subdirectory
			if not hasMultipleFiles and len(PATH_FILES) > 1:
				raise AssertionError('There are too many objects in path: '
					'{}'.format(PATH_DIR))

			# Make a copy of the current file in the parent directory
			elif len(PATH_FILES) == 1:
				for i in range(len(PATH_FILES)):
					IN_FILE = PATH_FILES[i]
					OUT_FILE = os.path.join(path, os.path.basename(PATH_FILES[i]))

					# Copy file to parent directory
					log_message('Writing to disk')
					
					with open(IN_FILE, 'rb') as f_in, \
							open(OUT_FILE, 'wb') as f_out:
						
						shutil.copyfileobj(f_in, f_out)
					
					# Log that back-up was successful
					log_message('Successfully backed up: {}'.format(IN_FILE))

				# return the file path if there is only one, else return the folder name
				return PATH_FILES[0] if not(hasMultipleFiles) else PATH_DIR

			# No files in current directory, so completely new run
			# OR
			# multiple files where there shouldn't be (input hasMultipleFiles == False, but count > 0)
			else:
				# get list of folder contents
				OLD_FILES = os.listdir(PATH_DIR)
				# filter out non-files
				OLD_FILES = list(filter(os.path.isfile, OLD_FILES))
				
				# if anything is left that wasn't caught above (count==1 or count > 1),
				# then raise error to stop run
				if len(OLD_FILES) > 0:
					raise AssertionError('We are missing our most current'
						'file in: {}'.format(PATH_DIR))
				else:
					# no files present here, so log message and return path to 'current' folder
					log_message('Initializing path: {} for organism: {}'.format(path, prefix))
					return PATH_DIR

		# Add files here in the future, make sure that the files follow the 
		# backup scheme in the above backup function
		files = {
			'tree': (os.path.join(RUNTIME_ARGS['dirpath'], 'tree'), False),
			'allele calls': (os.path.join(RUNTIME_ARGS['dirpath'], 'allele_calls'), True)
		}

		totalfiles = float(len(files))

		# For each file, back it up
		for (file, (path, hasMultipleFiles)) in files.items():
			log_message('Backing up if possible: {}'.format(path))
			
			# Set the runtime arguments
			RUNTIME_ARGS[file] = Backup(path, hasMultipleFiles)

			log_message('{} file backed up at {}'.format(file, RUNTIME_ARGS[file]))
			
		# Log successful back-up and return edited RUNTIME_ARGS
		log_message('Successfully validated source files')
		return RUNTIME_ARGS

	def ValidateLogFiles(RUNTIME_ARGS):
		"""
		ValidateLogFiles:  Add files paths to RUNTIME_ARGS dictionary for log files and Xcodes, if found
			Arguments:
				RUNTIME_ARGS:  dict --> original RUNTIME_ARGS dict declared above for editing and later return
			Returns:
				dict (edited RUNTIME_ARGS dict)
		"""
		# Create date-stamped file in change_log directory to record changed Allele Codes
		changeLogPath = '{}.tsv'.format(os.path.join(DATA_DIR, 
													'{}_nomenclature_logs'.format(prefix), 
													'change_log',
													Now().split(' ')[0]))
		if not os.path.exists(changeLogPath):
			print('Creating new change_log file ({})'.format(changeLogPath))
			with open(changeLogPath, 'w+') as f:
				f.write('{}\n'.format('\t'.join(['Key', 'OldValue', 'NewValue', 'ChangeType'])))
				f.close()
		# Add path to today's change log to to RUNTIME_ARGS
		RUNTIME_ARGS['changeLog'] = changeLogPath
		
		# Add path to Xcode file if found
		xCodeListPath = r'{}\{}\Xcodes.tsv'.format(RUNTIME_ARGS['logpath'], 'Xcodes')
		if os.path.exists(xCodeListPath):
			RUNTIME_ARGS['Xcodes'] = xCodeListPath
		
		# log successful validation and return edited RUNTIME_ARGS dictionary
		log_message('Successfully validated log files')
		return RUNTIME_ARGS
		
		
	# Verify source files and log directories exist
	ValidateDirectories()

	# Create a backup function to be run after the dialog box runs; 
	# returns edited RUNTIME_ARGS dictionary with file paths to use in algorithm
	RUNTIME_ARGS = ValidateSrcFiles(RUNTIME_ARGS)
	
	backup_function = partial(ValidateLogFiles, RUNTIME_ARGS)

	# return edited RUNTIME_ARGS dictionary with file paths to use in algorithm
	return backup_function


def Main(args):
	"""
	Main:  performs initial set-up, then runs algorithm using input arguments
	"""
	# Create data and logging directories
	initialize_srcfiles(DATA_DIR)
	initialize_logging(LOG_DIR)

	try:
		# Create a lock file to prevent concurrent data access from another user/instance
		lockfile_path = os.path.join(DATA_DIR,
									'{}_nomenclature_srcfiles'.format(prefix),
									'nomenclature.lock')

		if os.path.exists(lockfile_path):
			# file exists, so alert user and raise error to quit
			log_message('Lockfile was present on run')
				
			raise RuntimeError('Someone is running this already OR'
				' there was an error in the previous run')

		else:
			# lock file not present, so create it
			with open(lockfile_path, 'w') as f:
				pass

		# run Setup function to validate directories
		backupfnc = Setup()
		
		# get filled in RUNTIME_ARGS dictionary
		runtimeArgs = backupfnc()
			
		# Start the calculation
		Run(runtimeArgs)

	except:
		# log general error to file regardless of execution method
		log_exception('There was an error, please check the log file')
		

if __name__ == '__main__':
	
	# get inputs from user with defaults if argument(s) not provided
	inputArgs = argparse.ArgumentParser(prog=os.path.basename(__file__),
										formatter_class=argparse.RawDescriptionHelpFormatter,
										add_help=False, 
										description="Calculates nearest neighbor in hierarchical, single-linkage framework to assign that neighbor's code up to the corresponding distance threshold, or a new code if no matches found",
										epilog='Example usage:  python assignAlleleCodes_v3.6.py -a profiles.tsv -c locusNames.txt -d $(pwd) -p CAMP --nosave\n\n'
												'*Note:  locus order is determined by input config file and not saved to data directory.  The same file should be used '
												'in consecutive runs to ensure consistent locus-to-locus comparisons.')
	# required arguments:
	requiredArgs = inputArgs.add_argument_group('required arugments')
	
	# - tsv/csv of allele profiles (1st column is strain IDs, remaining are allele numbers with locus names as first line of file)
	requiredArgs.add_argument('-a', '--alleles', 
								nargs=1, 
								required=True, 
								type=str, 
								help='csv or tsv file.  Unique Keys/IDs in first column with PulseNet locus names as fields afterward, '
									'preceded with acceptable prefixes CAMP, EC, LMO, or SALM (e.g. EC_1, EC_2, EC_3, etc.')
	# - config file, containing names of core loci
	requiredArgs.add_argument('-c', '--config',
								nargs=1,
								required=True,
								type=str,
								help='text file containing names of core loci to match against --alleles file to denote which '
										'loci should be used for assignment')
	# - data directory holding existing tree, allele profiles, and logging directory
	requiredArgs.add_argument('-d', '--dataDir', 
								nargs=1, 
								required=True, 
								type=str,
								help='directory where logs and data files will be written if --nosave is not entered as input flag')
	# - prefix:  prefix to add to Allele Code with version if requested
	requiredArgs.add_argument('-p', '--prefix', 
								nargs=1, 
								required=True, 
								type=str, 
								choices=['CAMP', 'EC', 'LMO', 'SALM'],
								help='organism-specific prefix to be added to front of Allele Code')
	
	# optional arguments:
	optionalArgs = inputArgs.add_argument_group('optional arugments')
	optionalArgs.add_argument('-h', '--help', 
								action='help', 
								help='show this help message and exit')
	
	# - overwrite:  whether to create new tree file and overwrite historical allele calls or just export assigned codes to command-line
	optionalArgs.add_argument('--nosave', 
								action='store_true',
								help='whether results of run will overwrite historical allele profiles and nomenclature tree;  default = False if not provided')
	# - verbosity:  whether to output to only log file or also to terminal
	optionalArgs.add_argument('--verbose', 
								action='store_true',
								help="output what's written to log file if provided; default = False")
	# - output:  file to write results to for all new assignments and changed codes
	optionalArgs.add_argument('-o', '--output', 
								nargs=1, 
								type=str,
								help='output file (tsv or csv) into which results will be written (delimiter determined by input extension)')
	
	# get args from user into variables
	args = inputArgs.parse_args(sys.argv[1:])
	
	# check that required args are given
	if not args.alleles:
		print("New allele profiles (-a, --alleles) not provided")
		exit()
	if not args.config:
		print("Config file (--config) not provided.")
		exit()
	if not args.dataDir:
		print("Data directory (-d, --dataDir) not provided")
		exit();
	if not args.prefix:
		print("Organism prefix (-p, --prefix) required.  Options:  [CAMP, EC, LMO, SALM]")
		exit()
	
	# set global variables from user input
	newAllelesPath = args.alleles[0]
	if newAllelesPath.endswith('.csv'):
		delim = ','
	elif newAllelesPath.endswith('.tsv'):
		delim = '\t'
	else:
		print('Only csv and tsv files allowed for new allele profiles to be assigned an Allele Code.')
		exit()
		
	DATA_DIR = args.dataDir[0]
	prefix = args.prefix[0]
	configPath = args.config[0]
	with open(configPath, 'r') as c:
		coreLoci = [line.strip() for line in c if line.startswith(prefix)]
		
	# warn user if core locus names didn't load or prefix doesn't match
	if len(coreLoci)==0:
		print('Config file containing locus names not loaded.  Ensure organism prefix precedes locus names in config file and try again.')
		exit()
	
	
	# optional args
	if args.output:
		outputPath = args.output[0]
	nosave = args.nosave
	verbose = args.verbose
	
	# set log directory from input data directory
	LOG_DIR = os.path.join(DATA_DIR,
							'{}_nomenclature_logs'.format(prefix))
	
	# run Main function to set up and begin assignment
	Main({})
	
	
	