# remove VERTICES section from a VTK file 
# Forrest Sheng Bao http://fsbao.net
# GNU GPL v3.0 or later

def remove(inputvtk,outputvtk):
	"""Scan the inputvtk file line by line until the keyword VERTICES is found and remove corresponding lines
	"""
	fout=open(outputvtk, 'w')
	with open(inputvtk,'r') as fin:
		skip_lines = 0 
		for line in fin.readlines():
			if line[:8] == "VERTICES":
				line = line.split()
				skip_lines = 1 #int(line[1])  # the number of lines below should be skipped 
			else:
				if skip_lines > 0:
					skip_lines -= 1  # skip the line
				else:
					fout.write(line)

	fout.close()
	return 0

if __name__ == "__main__":
	import sys
	remove(sys.argv[1], sys.argv[2])
