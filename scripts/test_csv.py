import csv


pos_fil = open('/Users/elizabethdrenkard/external_data/CalCOFI/CalCOFIStaPosNDepth113.csv','rU')
rd_pos_fil = csv.reader(pos_fil, delimiter=',')  
mp=0
line_num = '93.3' 
for p_row in rd_pos_fil:
    if ((p_row[1] == line_num) and (p_row[2] == '28')):
       print p_row[3]
       print p_row[4]
    
