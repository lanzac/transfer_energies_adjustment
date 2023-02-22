#!/usr/bin/env python3

import re

# https://github.com/choderalab/ambermini/blob/master/share/amber/dat/leap/parm/parm99.dat

#*******************************Type**************************************#
file = open('nametotype.txt')
nametotypelines = file.readlines()
nametotype={}


pattern = re.compile('(\w+)(\t)(\w+)')

for nametotypeline in nametotypelines :
	if pattern.search(nametotypeline):
		m = pattern.match(nametotypeline)	    
		nametotype[m.group(1)]=m.group(3)
		
#print nametotype

file.close()

#*******************************Charge**************************************#
file = open('nametocharge.txt')
nametochargelines = file.readlines()
nametocharge={}


pattern = re.compile('(\w+)(\t)([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')

for nametochargeline in nametochargelines :
	if pattern.search(nametochargeline):
		m = pattern.match(nametochargeline)
		nametocharge[m.group(1)]=float(m.group(3))
		
#print nametocharge

file.close()
#*******************************Radius**************************************#
# Amber parm99 : https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html

file = open('typetoradius.txt')
typetoradiuslines = file.readlines()
typetoradius={}


pattern = re.compile('(\w+)(\t)([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')

for typetoradiusline in typetoradiuslines :
	if pattern.search(typetoradiusline):
		m = pattern.match(typetoradiusline)	    
		typetoradius[m.group(1)]=float(m.group(3))		
#print typetoradius

file.close()


#*******************************Mass**************************************#

file = open('typetomass.txt')
typetomasslines = file.readlines()
typetomass={}


pattern = re.compile('(\w+)(\t)([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')

for typetomassline in typetomasslines :
	if pattern.search(typetomassline):
		m = pattern.match(typetomassline)	    
		typetomass[m.group(1)]=float(m.group(3))		
#print typetomass

file.close()

#*******************************Epsilon**************************************#



file = open('typetoepsilon.txt')
typetoepsilonlines = file.readlines()
typetoepsilon={}


pattern = re.compile('(\w+)(\t)([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')

for typetoepsilonline in typetoepsilonlines :
	if pattern.search(typetoepsilonline):
		m = pattern.match(typetoepsilonline)	    
		typetoepsilon[m.group(1)]=float(m.group(3))		
#print typetoepsilon
file.close()
#*******************************Transfer**************************************#



file = open('typetotransfer.txt')
typetotransferlines = file.readlines()
typetotransfer={}

pattern = re.compile('(\w+)(\t)([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')

for typetotransferline in typetotransferlines :
	if pattern.search(typetotransferline):
		m = pattern.match(typetotransferline)	    
		typetotransfer[m.group(1)]=float(m.group(3))
		
#print typetotransfer

file.close()

#*****************************************************************************#
#*******************************Ajusted transfer**************************************#
file = open('nametotransfer_ajusted.txt')
nametotransferlines = file.readlines()
nametotransfer={}

pattern = re.compile('(\w+)(\t)([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')

for nametotransferline in nametotransferlines :
	if pattern.search(nametotransferline):
		m = pattern.match(nametotransferline)
		nametotransfer[m.group(1)]=float(m.group(3))
		
#print nametocharge

file.close()
#*******************************Radius**************************************#


file = open('amber_ajusted.ff', 'w')

file.write("#type\tcharge(e)\tradius(A)\tepsilon(kJ.mol-1)\tmass(Da)\ttransferIMP(kJ.mol-1.A-2\n")

list = nametotype.items()
list = sorted(list)
#print list
for item in list : 
	line = ""
	type = nametotype[item[0]]
	name = item[0]
	line += name
	 
	errorparsing = False
	if name in nametocharge: 
		#print str(nametocharge[name]) 
		line += "\t" + str(nametocharge[name])
	else : 
		print ("nametocharge has no key " + name )
		errorparsing = True

	if type in typetoradius: 
		#print str(typetoradius[type]) 
		line += "\t" + str(typetoradius[type])
	else : 
		print ("typetoradius has no key " + type )
		errorparsing = True

	if type in typetoepsilon: 
		#print str(typetoepsilon[type]) 
		line += "\t" + str(typetoepsilon[type])
	else : 
		print ("typetoepsilon has no key " + type )
		errorparsing = True

	if type in typetomass: 
		#print str(typetomass[type]) 
		line += "\t" + str(typetomass[type])
	else : 
		print ("typetomass has no key " + type )
		errorparsing = True

	# if type in typetotransfer: 
	# 	#print str(typetotransfer[type])
	# 	line += "\t" + str(typetotransfer[type])
	# else : 
	# 	print ("typetotransfer has no key " + type)
	# 	errorparsing = True

	if name in nametotransfer: 
		#print str(nametotransfer[name])
		line += "\t" + str(nametotransfer[name])
	else : 
		print ("nametotransfer has no key " + name)
		errorparsing = True


	if not errorparsing :  
		file.write(line+"\n")


file.close()


