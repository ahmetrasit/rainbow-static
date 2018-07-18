

def groupByIndex(list, index):
	types = {}
	for element in sorted(list, key=lambda x: x.split(':')[index]):
		type = element.split(':')[index]
		if type in types:
			types[type].append(element)
		else:
			types[type] = [element]
			
	return types
	
	
	

def parseFasta(file):
	f = open(file)
	exons = ''.join(f.readlines()).split('>')[1:]
	all = []
	for exon in exons:
		temp = exon.split('\n')
		header = ':'.join(temp[0].split('|'))
		seq = ''.join([x.strip() for x in temp[1:]])
		if '::' in header:
			pass
		else:
			all.append(header + ':' + seq)
	
	return all
	

			
def calculateThreshold(genomeSize, degreeLimit):
	return int( genomeSize * degreeLimit / (360*.98) )
	


def groupNonCodingElements(list, threshold):
	#first, group by type
	types = groupByIndex(list, 2)
	for type in types.keys():
		print type
		currType = types[type]
		#group by chromosome
		chrs = groupByIndex(currType, 1)
		for chr in chrs:
			print chr
			currChr = chrs[chr]
			#group by strand
			strands = groupByIndex(currChr, 3)
			for strand in strands:
				print strand
				curr = strands[strand]
				genes = findGenes(curr)
				
				file = open('./data/' + type + '_' + chr + '_' + strand + '.csv', 'w')
				output = 'start,end,amount,elements\n'
				
				regions, elementsInRegion = findFullRegions(genes, threshold)
				for i, region in enumerate(regions):
					if len(elementsInRegion[i]) > 0:
						output += ','.join([str(x) for x in region]) + ',' + str(len(elementsInRegion[i])) + ','  +';'.join(elementsInRegion[i]) + '\n'
				file.write(output)
				file.close()
				
				
				
				
				
def findGenes(list):
	genes = {}
	for exon in list:
		gene = exon.split(':')[0]
		if gene in genes:
			genes[gene].append(exon.replace(';', '_'))
		else:
			genes[gene] = [exon.replace(';', '_')]
	sortedGenes = []
	for gene in sorted(genes):
		sortedGene = sorted(genes[gene], key=lambda x: int(x.split(':')[5].split('_')[0]))
		start = int(sortedGene[0].split(':')[5].split('_')[0])
		end = int(sortedGene[-1].split(':')[6].split('_')[0])
		if (start > end):
			sortedGenes.append(str(end) + '|' + str(start) + '|' + '|'.join(sortedGene))
		else:
			sortedGenes.append(str(start) + '|' + str(end) + '|' + '|'.join(sortedGene))
			
	return sortedGenes
	
	
	
def findFullRegions(list, threshold):
	prev = 1
	regions = [[0,1]]
	elementsInRegion = [[]]
	for element in sorted(list, key=lambda x: int(x.split('|')[0]) ):
		start = int(element.split('|')[0])
		end = int(element.split('|')[1])
		if start - prev > threshold:
			regions.append([start,end])
			prev = end
			elementsInRegion.append([element])
		else:
			regions[-1][1] = end
			prev = end
			elementsInRegion[-1].append(element)
	
	return regions, elementsInRegion		
	
	
	
	
		
		
#chromosome IV:
all = parseFasta('../cel235ExonsBiomart.txt')
threshold =calculateThreshold(17493829, 0.1)
groupNonCodingElements(all, threshold)	
	
				