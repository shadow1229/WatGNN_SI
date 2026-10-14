import os,sys
trg_time = {}
with open('targets_cmp','r') as trg_file:
    trg_lines = trg_file.readlines()
    for trg_line in trg_lines:
        trg = trg_line.strip()
        trg_time[trg] = 0
        
trg_list = list(trg_time.keys())
trg_list.sort()

with open ('wkgb_010_time_freeze.log','r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('['):
            lsp = line.split()
            trg = lsp[1][:6]
            time = float(lsp[4].strip('s,'))
            if trg in trg_list:
                trg_time[trg] += time
                
out_f = open('superwater_time_summary.txt','w')
for trg in trg_list:
        out_f.write(f'{trg} {"%8.3f"%trg_time[trg]} \n')
