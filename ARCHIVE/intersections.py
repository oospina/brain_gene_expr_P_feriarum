# Python3 program to find common elements  
# in four lists using sets 

def IntersecOfSets(file1_array, file2_array, file3_array, file4_array): 
	# Converting the arrays into sets 
	s1 = set(file1_array) 
	s2 = set(file2_array) 
	s3 = set(file3_array) 
	s4 = set(file4_array)

	# Calculates intersections 
	set1_2 = s1.intersection(s2)
	set1_3 = s1.intersection(s3)
	set1_4 = s1.intersection(s4)

	set2_3 = s2.intersection(s3)    
	set2_4 = s2.intersection(s4)

	set3_4 = s3.intersection(s4)
 
	set1_2_set3 = set1_2.intersection(s3)
	set1_2_set4 = set1_2.intersection(s4)
	
	set1_3_set4 = set1_3.intersection(s4)
 
	set2_4_set3 = set2_4.intersection(s3)
 
	set1_2_set3_4 = set1_2.intersection(set3_4)
 
	# Converts resulting sets to lists and print 
	set1_2_list = list(set1_2) 
	print (file1, "vs", file2, ":")
	print (*set1_2_list, sep = "")

	set1_3_list = list(set1_3) 
	print (file1, "vs", file3, ":")
	print (*set1_3_list, sep = "")

	set1_4_list = list(set1_4) 
	print (file1, "vs", file4, ":")
	print (*set1_4_list, sep = "")

	set2_3_list = list(set2_3) 
	print (file2, "vs", file3, ":")
	print (*set2_3_list, sep = "")
 
	set2_4_list = list(set2_4) 
	print (file2, "vs", file4, ":")
	print (*set2_4_list, sep = "")
 
	set3_4_list = list(set3_4) 
	print (file3, "vs", file4, ":")
	print (*set3_4_list, sep = "")
 
	set1_2_set3_list = list(set1_2_set3) 
	print (file1,"-",file2, "vs", file3, ":")
	print (*set1_2_set3_list, sep = "")
 
	set1_2_set4_list = list(set1_2_set4) 
	print (file1,"-",file2, "vs", file4, ":")
	print (*set1_2_set4_list, sep = "")
 
	set1_3_set4_list = list(set1_3_set4) 
	print (file1,"-",file3, "vs", file4, ":")
	print (*set1_3_set4_list, sep = "")

	set2_4_set3_list = list(set2_4_set3) 
	print (file2,"-",file4, "vs", file3, ":")
	print (*set2_4_set3_list, sep = "")

	set1_2_set3_4_list = list(set1_2_set3_4) 
	print (file1,"-",file2, "vs", file3,"-",file4, ":")
	print (*set1_2_set3_4_list, sep = "")

	Allset_list = set1_2_list + set1_3_list + set1_4_list + set2_3_list + set2_4_list + set3_4_list + set1_2_set3_list + set1_2_set4_list + set1_3_set4_list + set2_4_set3_list + set1_2_set3_4_list

	s1unique = [item for item in s1 if item not in Allset_list]
	print (file1, "uniques:")
	print (*s1unique, sep = "")

	s2unique = [item for item in s2 if item not in Allset_list]
	print (file2, "uniques:")
	print (*s2unique, sep = "")

	s3unique = [item for item in s3 if item not in Allset_list]
	print (file3, "uniques:")
	print (*s3unique, sep = "")

	s4unique = [item for item in s4 if item not in Allset_list]
	print (file4, "uniques:")
	print (*s4unique, sep = "")

# Driver Code 
if __name__ == '__main__' : 

	file1 = 'topTags_FEMALES_SymVsAllo_ONLYIDs.txt'
	file2 = 'topTags_MALES_SymVsAllo_ONLYIDs.txt'
	file3 = 'topTags_SYMPATRY_FVsM_ONLYIDs.txt'
	file4 = 'topTags_ALLOPATRY_FVsM_ONLYIDs.txt'

	file1_array = []
	file2_array = []
	file3_array = []
	file4_array = []

with open(file1) as f1:
	for line1 in f1:
		file1_array.append(line1)

with open(file2) as f2:
	for line2 in f2:
		file2_array.append(line2)
  
with open(file3) as f3:
	for line3 in f3:
		file3_array.append(line3)
  
with open(file4) as f4:
	for line4 in f4:
		file4_array.append(line4)

# Calling Function 
IntersecOfSets(file1_array, file2_array, file3_array, file4_array)
