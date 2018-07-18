#WBGene00014450:tRNA:MtDNA:1:55:1:CAGTAAATAGTTTAATAAAAATATAGCATTTGGGTTGCTAAGATATTATTACTGA:MTCE.1


def findFullRegions(list, threshold):
	prev = 1
	regions = [[0,1]]
	elementsInRegion = [[]]
	for element in sorted(list, key=lambda x: int(x.split(':')[3]) ):
		start = int(element.split(':')[3])
		end = int(element.split(':')[4])
		if start - prev > threshold:
			regions.append([start,end])
			prev = end
			elementsInRegion.append([element])
		else:
			regions[-1][1] = end
			prev = end
			elementsInRegion[-1].append(element)
	
	return regions, elementsInRegion		
			
	
	
	
def groupByIndex(list, index):
	types = {}
	for element in sorted(list, key=lambda x: x.split(':')[index]):
		type = element.split(':')[index]
		if type in types:
			types[type].append(element)
		else:
			types[type] = [element]
			
	return types
	


def groupNonCodingElements(list, threshold):
	#first, group by type
	types = groupByIndex(list, 1)
	others = ['miRNA', 'tRNA', 'snRNA', 'lincRNA', 'antisense', 'snoRNA', 'rRNA']
	types['other'] = []
	for other in others:
		types['other'] += types[other]
		del types[other]
	print types.keys()	
	for type in types.keys():
		print type
		currType = types[type]
		#group by chromosome
		chrs = groupByIndex(currType, 2)
		for chr in chrs:
			print chr
			currChr = chrs[chr]
			#group by strand
			strands = groupByIndex(currChr, 5)
			for strand in strands:
				print strand
				curr = strands[strand]
				
				file = open('./data/' + type + '_' + chr + '_' + strand + '.csv', 'w')
				output = 'start,end,amount,elements\n'
				
				regions, elementsInRegion = findFullRegions(curr, threshold)
				for i, region in enumerate(regions):
					if len(elementsInRegion[i]) > 0:
						output += ','.join([str(x) for x in region]) + ',' + str(len(elementsInRegion[i])) + ','  +';'.join(elementsInRegion[i]) + '\n'
				file.write(output)
				file.close()
				
		
			
def calculateThreshold(genomeSize, degreeLimit):
	return int( genomeSize * degreeLimit / (360*.98) )
	
		
		
		
		
f = open('ncRNAs.parsed.txt')	
all = [x.strip() for x in f.readlines()]
threshold =calculateThreshold(17493829, 0.1)
groupNonCodingElements(all, threshold)	
	