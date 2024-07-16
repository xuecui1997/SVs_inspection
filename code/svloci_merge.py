#！svloci_merge.py
# -*- coding：utf-8 -*-
"""
Created on 2022/10/7 下午12:41

@author: cuixue

IDE: PyCharm 
"""
import sys
import getopt
import re
# from collections import Counter

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:t:o:",["infile=","outfile=","right="])
    except getopt.GetoptError:
        print('svloci_merge.py -i <posfile> -t <parameter> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('svloci_merge.py -i <posfile> -t <parameter> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inname = arg
        elif opt in ("-t", "--tfile"):
            tname = arg
        elif opt in ("-o", "--ofile"):
            oname = arg
    return inname,tname,oname
if __name__ == "__main__":
    inputname,typename,outputname = main(sys.argv[1:])
    inputfile = open(inputname)
    outputfile = open(outputname,'w')
    # usefile = open(usename,'w')
    # total = 0
    # overlap = 0
    # file1 = 0
    # file2 = 0
    # s = set()
    # l = []
    # # p = []
    # l_bef = ''
    # p_bef = 0
    s1 = ''
    pos1 = 0
    l1 = 0
    # a1 = ''
    old_start = 0
    new_start = 0
    old_chr = '0'
    line_num = 0
    type = set()
    svnum = 0
    # all_aligner = set()
    all_caller = set()
    last_end = 0
    tmp_end = set()
    st_set = []
    while 1:
        lines = inputfile.readlines(10000)
        if not lines:
            break
        for line1 in lines:
            line1 = line1.rstrip()
            cut1 = line1.strip().split('\t')
            if re.search('^#',line1):
                outputfile.write(line1 + '\n')
            else:
                now_chr = cut1[0]
                # print(cut1)
                now_start = cut1[1]
                now_end = cut1[2]
                svtype = cut1[3]
                caller = cut1[4]
                line_num += 1
                # if 'INS' in svtype:
                #     if line_num == 1:
                #         old_chr = now_chr
                #         old_start = now_start
                #         type.add(svtype)
                #         # all_aligner.add(aligner)
                #         all_caller.add(caller)
                #     else:
                #         if int(now_start) - int(old_start) > int(typename):
                #             type_end = ''
                #             aligner_end = ''
                #             caller_end = ''
                #             n = 0
                #             for i in type:
                #                 n += 1
                #                 if n == 1:
                #                     type_end = i
                #                 else:
                #                     type_end = type_end +'+'+ i
                #             # for j in all_aligner:
                #             #     aligner_end = all_aligner + '-' + j
                #             for k in all_caller:
                #                 n += 1
                #                 if n ==1:
                #                     caller_end = k
                #                 else:
                #                     caller_end = caller_end + '+' + k
                #                 outputfile.write(old_chr + '\t' + end_start + '\t' + old_end + '\t' + type_end + '\t' +
                #                                  caller_end + '\n')
                #             # outputfile.write(old_chr+'\t'+old_start+'\t'+end+'\t'+type_end+'\t'+aligner_end+'\t' + caller_end + '\n')
                #
                #             old_chr = now_chr
                #             old_start = now_start
                #             end = now_end
                #             type = set()
                #             # all_aligner = set()
                #             all_caller = set()
                #             type.add(svtype)
                #             # all_aligner.add(aligner)
                #             all_caller.add(caller)
                #         else:
                #             if old_chr == now_chr:
                #
                #                 old_chr = now_chr
                #             else:
                #                 print(line1)
                #             end_start = old_start
                #             old_start = now_start
                #             end = now_end
                #             type.add(svtype)
                #             # all_aligner.add(aligner)
                #             all_caller.add(caller)
##############other svtype merge maotheds:
#1       167747  172630  INS01   assemblytics
# 1       167752  172624  INS11   pbsv.pbmm2
# 1       167752  172624  none    pbsv.ngmlr
# 1       167752  172625  INS11   syri
# 1       167758  172619  INS11   svmu
#2       395883  395947  DEL11   svision.pbmm2
#2	22199293	22199364	none	pbsv.pbmm2
                if line_num == 1:
                    old_chr = now_chr
                    old_start = now_start
                    end_start = now_start
                    old_start_1 = now_start
                    old_end = now_end
                    old_type = svtype
                    old_caller = caller
                    old_end_1 = now_end
                    s_t = now_start +':' + svtype
                    st_set.append(s_t)
                    tmp_end.add(old_end_1)
                    type.add(svtype)
                    all_caller.add(caller)
                else:
# 1	10708170	10708248	INS01
# 1	10708186(old-s)	10708264(old-e)	INS01 input

# 1	10709557(now-s--olds)	10719220(now-e---olde)	DEL01
# 1	10709557(nows)	10719220(nowe)	none
# 1	10709557	10719220	none
# 1	10709558	10719219	none
# 1	10709559	10719221	none
# 1	10709560	10720615	none
# 1	10711293	10711394	none
# 1	10711418(olds)	10711531(olde)	DEL01

# 1	10712315(nows)	10978066(nowe)	INS01+DEL01+DUP+INV input

# 1	10719683	10719762	DEL01

                    if int(now_end) - int(now_start) > 200000 and \
                            int(now_start) - int(old_start) > 0:
                            # print(line1)
                        if int(now_end) - int(old_end) > 200000:
                            # print(line1)
                            svnum += 1
                            outputfile.write(line1+'\n')
                            ##test:
                            # print(line1)
                            # print(old_chr + '\t' + end_start + '\t' + old_end + '\t' + type_end + '\t' + caller_end + '\n')
                            juge = set()
                            for i in tmp_end:
                                if int(now_end) - int(i) > 0:
                                    juge.add('+')
                                else:
                                    juge.add('-')
                            if '-' not in juge:
                                type_end = ''
                                # aligner_end = ''
                                caller_end = ''
                                n = 0
                                for i in type:
                                    n += 1
                                    if n == 1:
                                        type_end = i
                                    else:
                                        type_end = type_end + '-' + i
                                # for j in all_aligner:
                                #     aligner_end = all_aligner + '-' + j
                                m = 0
                                for k in all_caller:
                                    m += 1
                                    if m == 1:
                                        caller_end = k
                                    else:
                                        caller_end = caller_end + '-' + k
                            svnum += 1
                            old_end = max(tmp_end)
                            outputfile.write(old_chr + '\t' + end_start + '\t' + old_end + '\t' + type_end + '\t' + caller_end + '\n')
                            old_start = 0
                            new_start = 0
                            old_chr = '0'
                            line_num = 0
                            type = set()
                            # svnum = 0
                            # all_aligner = set()
                            all_caller = set()
                            last_end = 0
                            tmp_end = set()
                            st_set = []

                        else:
                            continue
                    else:
                        if abs(int(now_start) - int(old_start)) > int(typename) and \
                                abs(int(now_end) - int(old_end)) > int(typename):
                            #
                            #### same chr:
                            if int(now_start) - int(old_start) > 0:
                                # print(tmp_end)
                                juge = set()
                                for i in tmp_end:
                                    if int(now_end) - int(i) >0:
                                        juge.add('+')
                                    else:
                                        juge.add('-')
                                if '-' not in juge:
                                    type_end = ''
                                    # aligner_end = ''
                                    caller_end = ''
                                    n = 0
                                    for i in type:
                                        n += 1
                                        if n == 1:
                                            type_end = i
                                        else:
                                            type_end = type_end + '-' + i
                                    # for j in all_aligner:
                                    #     aligner_end = all_aligner + '-' + j
                                    m = 0
                                    for k in all_caller:
                                        m += 1
                                        if m == 1:
                                            caller_end = k
                                        else:
                                            caller_end = caller_end + '-' + k
                                    svnum += 1
                                    ####juge true start and end:

                                    # p = st_set[0].split(':')
                                    # p_e = st_set[len(st_set)-1].split(':')
                                    # if p[1] == 'none' and len(type) >= 2:
                                    #     cir = 0
                                    #     for i in st_set:
                                    #         cir += 1
                                    #         i = st_set[cir].split(':')
                                    #         if i[1] == 'none':
                                    #             continue
                                    #         else:
                                    #             p1 = st_set[cir].split(':')
                                    #             end_start = p1[0]
                                    #             break

                                    # if p[1] == 'none' and len(type) >= 2:
                                    #     p1 = st_set[1].split(':')
                                    #     end_start = p1[0]
                                    #     # print(end_start)
                                    # if p_e[1] == 'none' and len(type) >= 2:
                                    #     p_e1 = st_set[len(st_set) - 2].split(':')
                                    #     old_end = p_e1[0]
                                    ###### write out the true svsite:
                                    old_end = max(tmp_end)
                                    outputfile.write(old_chr + '\t' + end_start + '\t' + old_end + '\t' + type_end + '\t' + caller_end + '\n')
                                    old_chr = now_chr
                                    old_start = now_start
                                    old_end = now_end
                                    old_start_1 = now_start
                                    end_start = now_start
                                    old_end_1 = now_end
                                    old_type = svtype
                                    old_caller = caller
                                    tmp_end = set()
                                    tmp_end.add(old_end_1)
                                    type = set()
                                    all_caller = set()
                                    st_set =[]
                                    s_t = now_start + ':' + svtype
                                    st_set.append(s_t)
                                    type.add(svtype)
                                    all_caller.add(caller)
                                else:
                                    juge1 = set()
                                    for i in type:
                                        if 'INV' in i:
                                            juge1.add('break')
                                        else:
                                            juge1.add('continue')
                                    if 'break' in juge1:
                                        p = st_set[0].split(':')
                                        p_e = st_set[len(st_set) - 1].split(':')

                                        # if p[1] == 'none' and len(type) >= 2:
                                        #     p1 = st_set[1].split(':')
                                        #     end_start = p1[0]
                                        #     # print(end_start)
                                        # if p_e[1] == 'none' and len(type) >= 2:
                                        #     p_e1 = st_set[len(st_set) - 2].split(':')
                                        #     old_end = p_e1[0]
                                        svnum += 1
                                        old_end = max(tmp_end)
                                        svinfo = old_chr+'\t'+old_start+'\t'+old_end+'\t'+old_type+'\t'+old_caller+'\n'
                                        outputfile.write(svinfo)

                                        old_chr = now_chr
                                        old_start = now_start
                                        old_end = now_end
                                        old_start_1 = now_start
                                        end_start = now_start
                                        old_end_1 = now_end
                                        old_type = svtype
                                        old_caller = caller
                                        tmp_end = set()
                                        tmp_end.add(old_end_1)
                                        type = set()
                                        # all_aligner = set()
                                        all_caller = set()
                                        st_set = []
                                        s_t = now_start + ':' + svtype
                                        st_set.append(s_t)
                                        type.add(svtype)
                                        # all_aligner.add(aligner)
                                        all_caller.add(caller)
                                    else:
                                        old_chr = now_chr
                                        # else:
                                        #     print(line1)
                                        end_start = old_start_1
                                        old_start = now_start
                                        old_end = now_end
                                        old_type = svtype
                                        old_caller = caller
                                        tmp_end.add(old_end)
                                        type.add(svtype)
                                        s_t = now_start + ':' + svtype
                                        st_set.append(s_t)
                                        # all_aligner.add(aligner)
                                        all_caller.add(caller)

                            else:
                                type_end = ''
                                # aligner_end = ''
                                caller_end = ''
                                n = 0
                                for i in type:
                                    n += 1
                                    if n == 1:
                                        type_end = i
                                    else:
                                        type_end = type_end + '-' + i
                                # for j in all_aligner:
                                #     aligner_end = all_aligner + '-' + j
                                m = 0
                                for k in all_caller:
                                    m += 1
                                    if m == 1:
                                        caller_end = k
                                    else:
                                        caller_end = caller_end + '-' + k
                                svnum += 1
                                p = st_set[0].split(':')
                                p_e = st_set[len(st_set) - 1].split(':')

                                # if p[1] == 'none' and len(type) >= 2:
                                #     p1 = st_set[1].split(':')
                                #     end_start = p1[0]
                                #     # print(end_start)
                                # if p_e[1] == 'none' and len(type) >= 2:
                                #     p_e1 = st_set[len(st_set) - 2].split(':')
                                #     old_end = p_e1[0]
                                old_end = max(tmp_end)
                                outputfile.write(old_chr + '\t' + end_start + '\t' + old_end + '\t' + type_end + '\t' + caller_end + '\n')
                                old_chr = now_chr
                                old_start = now_start
                                old_end = now_end
                                old_start_1 = now_start
                                end_start = now_start
                                old_end_1 = now_end
                                old_type = svtype
                                old_caller = caller
                                tmp_end = set()
                                tmp_end.add(old_end_1)
                                type = set()
                                # all_aligner = set()
                                all_caller = set()
                                st_set = []
                                s_t = now_start + ':' + svtype
                                st_set.append(s_t)
                                type.add(svtype)
                                # all_aligner.add(aligner)
                                all_caller.add(caller)
                        else:
                            # if abs(int(now_start) - int(old_start)) < int(typename):
                            #     if int(now_end) - int(old_end) < 0:
                            #         print(line1)
                            # if now_chr == old_chr:

                            old_chr = now_chr
                            # else:
                            #     print(line1)
                            end_start = old_start_1
                            old_start = now_start
                            old_end = now_end
                            old_type = svtype
                            old_caller = caller
                            s_t = now_start + ':' + svtype
                            st_set.append(s_t)
                            tmp_end.add(old_end)
                            type.add(svtype)
                            # all_aligner.add(aligner)
                            all_caller.add(caller)






print('sv-merge_num:',svnum)
inputfile.close()
outputfile.close()




                # s.add(start)