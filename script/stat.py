#!/usr/bin/python
gammalist = ['0.1','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0']#,'4.0','10','20','40','50','100','200','300','1000','10000']
tlist = ['0.1','0.3']#,'0.05']
limit = 0.003
tmp = 0.123
for i in gammalist:
	for j in tlist:
		filename = 'corfile/'+j+'-'+i+'-mdcor'
		f = open(filename)
		for line in f.readlines():
			linelist = line.split()
			if abs(tmp-float(linelist[2]))<limit:
				s = ','.join(item for item in linelist)
				s = j + ',' + i + ',' + s + '\n'
				with open(j+'.csv','ab+') as ff:
					ff.write(s)
#				print s
				break
			tmp = float(linelist[2])
		
