#！SV_custom.py
# -*- coding：utf-8 -*-
"""
Created on 8/29/23 11:28 PM

@author: cuixue

IDE: PyCharm 
"""
import re
import argparse


def main(argparse):
    parser = argparse.ArgumentParser(description='SV detection')

    parser.add_argument('--input', '-i', help='inputfile', \
                        type=str)
    parser.add_argument('--expand', '-e', help='output: 左右扩展500bp 用作 INS 和 DEL 的侧翼序列定位，同时调整 INS 的起始位点==终止位点', \
                        type=str)
    parser.add_argument('--filter', '-f', help='output: 过滤后的位点文件', \
                        type=str)
    parser.add_argument('--exnum', '-n', help='para: 左右扩展的 bp 数量，eg：500', \
                        type=int)
    parser.add_argument('--bedfile', '-b', help='input: ref 位点提取 fq 重新比对到 qry.fa 上的全部的 bedfile 位置', \
                        type= str)
    parser.add_argument('--qrypos', '-q', help='output: 输出对应的 qry 上的 SV 的位置', \
                        type=str)
    parser.add_argument('--refchr', '-rc', help='para: 参考基因组的染色体', \
                        type=str)
    parser.add_argument('--refs', '-rs', help='para: 参考基因组的对应的 SV 的起始位点', \
                        type=str)
    parser.add_argument('--refe', '-re', help='para: 参考基因组的对应的 SV 的终止位点', \
                        type=str)
    parser.add_argument('--refqry', '-p', help='input: bedtools 得到的 ref 对应 qry 的位置', \
                        type=str)
    parser.add_argument('--mumex', '-x', help='output: 给到 mummer 扩充比对的位置', \
                        type=str)
    parser.add_argument('--unex', '-u', help='para: 判断比对到 qry 上 qry 长度终止位点不合理的阈值', \
                        type=int)
    parser.add_argument('--reffolder', '-rf', help='para: mummer批量的所有文件夹的名字', \
                        type=str)
    parser.add_argument('--linenumber', '-ln', help='para: mummer批量的所有对应的结果 mcroods 需要的行数', \
                        type=str)
    parser.add_argument('--mcoords', '-mc', help='input: mummer批量的所有对应的结果 mcroods文件', \
                        type=str)
    parser.add_argument('--indelfile', '-in', help='output: 得到了 ins 和 del 的对应的 type1+type2 的标准表 pre', \
                        type=str)
    parser.add_argument('--manual', '-ma', help='output: 需要 check 的连在一起的 type2+TR 等有 overlap 的 mcoods 文件和 sv 位点', \
                        type=str)

    parser.add_argument('--model', '-m', help='para: 执行哪个模块', \
                        type=str)
    ####
    parser.add_argument('--trfile', '-tr', help='output:  能够判断出可以输出一段范围的 TR', \
                        type=str)
    ###
    parser.add_argument('--invfile', '-inv', help='output:  能够判断出可以输出一段范围的 TR', \
                        type=str)
    ####
    parser.add_argument('--stand', '-sd', help='output:  输出为标准表格的文件', \
                        type=str)
    ####
    parser.add_argument('--vcffile', '-vcf', help='output:  输出为标准表格转后的 vcf 统一格式文件', \
                        type=str)
    #correctname
    parser.add_argument('--correctname', '-co', help='output:  输出为标准表格 更正一些位点后的最终标准表格', \
                        type=str)
    ###常规 output
    parser.add_argument('--output', '-out', help='output:  输出为常规的 output', \
                        type=str)
    args = parser.parse_args()
    return args.input, args.expand, args.filter, args.exnum, args.bedfile, args.qrypos, args.refchr, args.refs, \
           args.refe, args.refqry, args.mumex, args.unex, args.reffolder, args.linenumber, args.mcoords, \
           args.indelfile, args.manual,\
           args.model, args.trfile, args.invfile, args.stand, args.vcffile, args.correctname, args.output

if __name__ == '__main__':
    input, expandname, filtername, numname, bedname, qryname, refchrget, refsget, refeget, posname, mumname, unexpect,\
    refname,linename,mcoordname,insdelname,manualname,model,trname,invname,standname,vcfname,correctname,outputname \
        = main(argparse)
    # inputfile = open(input)

    if model == 'expand':
        inputfile = open(input)
        exfile = open(expandname, 'w')
        while True:
            lines = inputfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split(' ')
                    # print(cut)
                    refchr = cut[0]
                    refs = cut[1]
                    refe = cut[2]
                    type = cut[3]
                    aligner = cut[4]
                    if 'INS' in type and 'DEL' not in type:
                        finals = int(refs) - numname + 1
                        finale = int(refs) + numname
                    else:
                        finals = int(refs) - numname + 1
                        finale = int(refe) + numname - 1
                    length = finale-finals+1
                    exfile.write(refchr+'\t'+str(finals)+'\t'+str(finale)+'\t'+type+'\t'+aligner+'\t'+str(length)+'\n')
        inputfile.close()
        exfile.close()

    if model == 'bedtopos':
        bedfile = open(bedname)
        qryfile = open(qryname, 'a')

        l = 0
        qryall_s = []
        qryall_e = []
        qrychrall=set()
        qrychrchange ={'Chr01':'1','Chr02':'2','Chr03':'3','Chr04':'4','Chr05':'5','Chr06':'6','Chr07':'7',\
                       'Chr08':'8','Chr09':'9','Chr10':'10','Chr11':'11','Chr12':'12'}
        while True:
            lines1 = bedfile.readlines(10000)
            if not lines1:
                break
            else:
                for line1 in lines1:
                    line1 = line1.rstrip()
                    cut1 = line1.strip().split('\t')
                    l += 1
                    tmpqrychr = cut1[0]
                    tmpqrys = int(cut1[1])
                    tmpqrye = int(cut1[2])
                    qrychr = qrychrchange[tmpqrychr]
                    if qrychr in qrychrall:
                        qryall_s.append(tmpqrys)
                        qryall_e.append(tmpqrye)
                    else:
                        if l == 1:
                            if refchrget == qrychr:
                                finalqrychr = tmpqrychr
                                qrychrall.add(qrychr)
                                qryall_s.append(tmpqrys)
                                qryall_e.append(tmpqrye)
                            else:
                                continue
                        else:
                            if refchrget == qrychr:
                                # print(refchrget)
                                finalqrychr = tmpqrychr
                                qrychrall.add(qrychr)
                                qryall_s.append(tmpqrys)
                                qryall_e.append(tmpqrye)
                            else:
                                continue
        final_qrys = min(qryall_s)
        final_qrye = max(qryall_e)
        final = [refchrget,refsget,refeget,finalqrychr,str(final_qrys),str(final_qrye)]
        qryfile.write('\t'.join(final)+'\n')

        bedfile.close()
        qryfile.close()

    if model == 'qrychange':
        posfile = open(posname)
        mumfile = open(mumname, 'w')
        qrysall = []
        while True:
            lines = posfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    refchr = cut[0]
                    refs = cut[1]
                    refe = cut[2]
                    qrychr = cut[3]
                    qrys = cut[4]
                    qrye = cut[5]

                    if qrys not in qrysall:
                        qrysall.append(qrys)
                        reflength = int(refe) - int(refs) + 1
                        qrylength = int(qrye) - int(qrys) + 1
                        # if reflength == 1367450:
                        #     print(line)
                        if int(qrye) - int(qrys) > unexpect:
                            # mumfile.write(line+'\t'+str(int(qrye) - int(qrys)+1)+'\n')
                            qrye = int(qrys)+reflength
                            allex = ['0','0','0','0']
                        else:
                            if qrylength > reflength:
                                diff = qrylength-reflength
                                juge = 'qry'
                            else:
                                diff = reflength - qrylength
                                juge = 'ref'
                            difflen = len(str(diff))
                            diff_one = str(diff)[0]
                            n = difflen-1
                            exnum = int(int((diff_one+'0'*n))/2)
                            if juge == 'qry':
                                #refs- refe+ qrys- qrye+
                                allex = [str(exnum),str(exnum),'0','0']
                            else:
                                allex = ['0','0',str(exnum),str(exnum)]
                        exblock = '\t'.join(allex)
                        final = refchr+'\t'+refs+'\t'+refe+'\t'+\
                            qrychr+'\t'+qrys+'\t'+str(qrye)+'\t'+\
                            exblock+'\n'
                        mumfile.write(final)
                    else:
                        continue
        posfile.close()
        mumfile.close()

    if model == 'mcoordsw':
        #input: refchr_refs_refe line_num(1,2,4,5,6)
        #outout: stander_excel
        line_num = linename.split(',')
        numall = []
        for n in line_num:
            numall.append(int(n))
        def output(list,  overjuge,nowjuge, linenum,outfile, refname,numall,manualfile,lastsvtype,trfile):
            try:
                refchr = list[0]
                refsv_S = list[1]
                refsv_E = list[2]
                qrychr = list[3]
                qrysv_S = list[4]
                qrysv_E = list[5]
                svtype = list[6]
                overlap = list[7]
                overwho = list[8]
                sub = list[9]
                refori_s = list[10]
                f_ref_s1 = list[11]
                f_ref_e1 = list[12]
                f_ref_s2 = list[13]
                f_ref_e2 = list[14]
                qryori_s = list[15]
                f_qry_s1 = list[16]
                f_qry_e1 = list[17]
                f_qry_s2 = list[18]
                f_qry_e2 = list[19]
            except IndexError:
                pass
            # print(linenum)
            # print(overjuge,nowjuge)
            if linenum == 3:
                if overjuge == 'false' and nowjuge == 'true':
                    f_ref_s2 = '0'
                    f_ref_e2 = '0'
                    f_qry_s2 = '0'
                    f_qry_e2 = '0'
                    final = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                             svtype, overlap, overwho, sub, \
                             refori_s, \
                             f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                             qryori_s, \
                             f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                    f = '\t'.join(final)
                    outfile.write(f + '\n')
                if overjuge == 'true' and nowjuge == 'false':
                    if lastsvtype == 'TR':
                        final = '\t'.join(list)
                        trfile.write(final+'\n')
                        # f = [refname, fn, svtype, refori_s, refbs, refbe, qryori_s, qrybs, qrybe]
                        # manualfile.write('\t'.join(f) + '\n')
                    else:
                        final = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                 svtype, overlap, overwho, sub, \
                                 refori_s, \
                                 f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                 qryori_s, \
                                 f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                        f = '\t'.join(final)
                        outfile.write(f + '\n')
                if overjuge == 'true' and nowjuge == 'true':
                    if lastsvtype == 'TR':
                        co = 'TR'
                    else:
                        co = 'Type2'
                    # print(linenum)
                    f = refname + '\t' +str(numall[linenum - 2 - 1]) +','+str(numall[linenum - 1 - 1])+'\t' + co + \
                        '\t' + '\n'
                    manualfile.write(f)
                if overjuge == 'false' and nowjuge == 'false':
                    final = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                             svtype, overlap, overwho, sub, \
                             refori_s, \
                             f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                             qryori_s, \
                             f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                    f = '\t'.join(final)
                    outfile.write(f + '\n')
            else:
                if overjuge == 'false' and nowjuge == 'true':
                    f_ref_s2 = '0'
                    f_ref_e2 = '0'
                    f_qry_s2 = '0'
                    f_qry_e2 = '0'
                    f_ref_s1 = '0'
                    f_ref_e1 = '0'
                    f_qry_s1 = '0'
                    f_qry_e1 = '0'
                    final = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                             svtype, overlap, overwho, sub, \
                             refori_s, \
                             f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                             qryori_s, \
                             f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                    f = '\t'.join(final)
                    outfile.write(f + '\n')
                if overjuge == 'true' and nowjuge == 'false':
                    if lastsvtype == 'TR':
                        final = '\t'.join(list)
                        trfile.write(final+'\n')
                    else:
                        final = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                 svtype, overlap, overwho, sub, \
                                 refori_s, \
                                 f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                 qryori_s, \
                                 f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                        f = '\t'.join(final)
                        outfile.write(f + '\n')
                if overjuge == 'true' and nowjuge == 'true':
                    # n1 = numall[linenum]
                    # n2 = numall[linenum - 1 -1]
                    # manualfile.write(refname + '\t' + str(n1) + ',' + str(n2) +'\t'+'type2'+ '\n')
                    if lastsvtype == 'TR':
                        co = 'TR'
                    else:
                        co = 'Type2'
                    f = refname + '\t' +str(numall[linenum - 2 - 1]) +','+str(numall[linenum - 1 - 1])+'\t' + co + \
                        '\t' + '\n'
                    manualfile.write(f)
                if overjuge == 'false' and nowjuge == 'false':
                    f_ref_s1 = '0'
                    f_ref_e1 = '0'
                    f_qry_s1 = '0'
                    f_qry_e1 = '0'
                    final = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                             svtype, overlap, overwho, sub, \
                             refori_s, \
                             f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                             qryori_s, \
                             f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                    f = '\t'.join(final)
                    outfile.write(f + '\n')

        mcoordfile = open(mcoordname)
        manualfile = open(manualname, 'a')
        insdelfile = open(insdelname,'a')
        trfile = open(trname,'a')
        invfile = open(invname, 'a')
        l = 0
        real = 0
        truestore = []
        lastsvtype = ''
        overfind = 'no'
        while True:
            lines = mcoordfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    qry_opps = cut[0]
                    qry_oppe = cut[1]
                    ref_opps = cut[2]
                    ref_oppe = cut[3]
                    # print(line)
                    refinfo = cut[12].split(':')
                    refchr = refinfo[0]
                    refori_s = refinfo[1].split('-')[0]
                    qryinfo = cut[11].split(':')
                    qrychr = qryinfo[0]
                    qryori_s = qryinfo[1].split('-')[0]
                    l += 1

                    if l in numall:
                        if int(ref_oppe) < int(ref_opps):
                            invfile.write(refname + '\n')
                            break
                        real += 1
                        # print(numall)
                        if numall.index(l) == 0:

                            ref1 = int(ref_oppe)
                            qry1 = int(qry_oppe)
                            f_ref_s1 = ref_opps
                            f_qry_s1 = qry_opps
                            copy = 1
                        else:
                            ref2 = int(ref_opps)
                            qry2 = int(qry_opps)
                            if ref2 < ref1:#over-type1+2-INS
                                overlen = abs(ref2 - ref1) + 1
                                if qry1 < qry2:#type1
                                    # print(line)
                                    # print(numall.index(real))
                                    if overlen < 45:
                                        # print('1')
                                        svtype = 'INS'
                                        overlap = '0'
                                        overwho = '0'
                                        sub = '0'
                                        refsv_S = str(ref1)
                                        refsv_E = str(ref1)
                                        qrysv_S = str(qry1+1)
                                        qrysv_E = str(qry2+overlen -1-1)
                                        f_ref_e1 = refsv_S
                                        f_qry_e1 = str(qry1)
                                        f_ref_s2 = refsv_S
                                        f_ref_e2 = ref_oppe
                                        f_qry_s2 = str(qry2+overlen-1)
                                        f_qry_e2 = qry_oppe
                                        overfind = 'yes'
                                        newqry1 = str(int(qrysv_E)+1)
                                        newref1 = refsv_E
                                        if real >= 3:
                                            # print(line)
                                            # print(tmp_record)
                                            if truestore[-2:] == ['true', 'true']:
                                                if lastsvtype == 'TR':
                                                    co = 'TR'
                                                else:
                                                    co = 'Type2'
                                                f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                                manualfile.write(f)
                                            else:
                                                output(tmp_record, overjuge, 'false', real, insdelfile,
                                                   refname, numall, manualfile,lastsvtype,trfile)
                                            if real == len(numall):
                                                finaljuge = 'false'
                                                truestore.append(finaljuge)
                                            else:
                                                overjuge = 'false'
                                                truestore.append(overjuge)
                                        else:

                                            overjuge = 'false'
                                            truestore.append(overjuge)

                                        tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                      svtype, overlap, overwho, sub, \
                                                      refori_s, \
                                                      f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                      qryori_s, \
                                                      f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                        lastsvtype = 'INS'
                                    else:#type2:
                                        svtype = 'INS'
                                        overfind = 'no'
                                        overlap = '1'
                                        overwho = 'qry'
                                        sub = '0'
                                        refsv_S = str(ref2)
                                        refsv_E = str(ref1)
                                        qrysv_S = str(qry1)
                                        qrysv_E = str(qry2)
                                        f_ref_e1 = 'uncertain'
                                        f_qry_e1 = 'uncertain'
                                        f_ref_s2 = 'uncertain'
                                        f_ref_e2 = ref_oppe
                                        f_qry_s2 = 'uncertain'
                                        f_qry_e2 = qry_oppe
                                        if real >= 3:
                                            if truestore[-2:] == ['true', 'true']:
                                                if lastsvtype == 'TR':
                                                    co = 'TR'
                                                else:
                                                    co = 'Type2'
                                                f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                                manualfile.write(f)
                                            else:
                                                output(tmp_record, overjuge, 'true', real, insdelfile,
                                                       refname, numall, manualfile,lastsvtype,trfile)
                                            if real == len(numall):
                                                finaljuge = 'true'
                                                truestore.append(finaljuge)
                                            else:
                                                overjuge = 'true'
                                                truestore.append(overjuge)
                                        else:

                                            overjuge = 'true'
                                            truestore.append(overjuge)

                                        tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                      svtype, overlap, overwho, sub, \
                                                      refori_s, \
                                                      f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                      qryori_s, \
                                                      f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                        lastsvtype = 'INS'

                                else:
                                    # print(real)
                                    svtype = 'TR'
                                    overfind = 'no'
                                    # n1 = numall[real - 1]
                                    # n2 = numall[real - 1 - 1]
                                    refbs = f_ref_s1
                                    refbe = ref_oppe
                                    qrybs = f_qry_s1
                                    qrybe = qry_oppe
                                    reflenth1 = ref1 - ref2 +1
                                    qrylength1 = qry1 - qry2 + 1
                                    if reflenth1 > qrylength1:
                                        trtype = 'INS'
                                    else:
                                        trtype = 'DEL'
                                    # fn = str(n2) + ',' + str(n1)
                                    if real >= 3:
                                        if truestore[-2:] == ['true', 'true']:
                                            if lastsvtype == 'TR':
                                                co = 'TR'
                                            else:
                                                co = 'Type2'
                                            f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                            manualfile.write(f)
                                        else:
                                            output(tmp_record, overjuge, 'true', real, insdelfile,
                                                   refname, numall, manualfile, lastsvtype,trfile)
                                        if real == len(numall):
                                            finaljuge = 'true'
                                            truestore.append(finaljuge)
                                        else:
                                            overjuge = 'true'
                                            truestore.append(overjuge)
                                    else:
                                        overjuge = 'true'
                                        truestore.append(overjuge)
                                    # tmp_record = [refname,refchr,qrychr,trtype,refori_s,refbs,refbe,qryori_s,qrybs,
                                    #               qrybe]
                                    tmp_record = [refchr, str(ref2), str(ref1),'1', qrychr, str(qry2), str(qry1), \
                                                  '1',trtype, refori_s, \
                                                 refbs, 'uncertain', 'uncertain', refbe, \
                                                  qryori_s, \
                                                  qrybs,'uncertain', 'uncertain',qrybe]
                                    lastsvtype = 'TR'

                            elif ref2 == ref1:
                                #type1 =INS
                                if qry1 < qry2:
                                    svtype = 'INS'
                                    overfind = 'no'
                                    overlap = '0'
                                    overwho = '0'
                                    sub = '0'
                                    refsv_S = str(ref1)
                                    refsv_E = str(ref1)
                                    qrysv_S = str(qry1 + 1)
                                    qrysv_E = str(qry2 - 1)
                                    f_ref_e1 = refsv_S
                                    f_qry_e1 = str(qry1)
                                    f_ref_s2 = refsv_S
                                    f_ref_e2 = ref_oppe
                                    f_qry_s2 = str(qry2)
                                    f_qry_e2 = qry_oppe
                                    if real >= 3:
                                        if truestore[-2:] == ['true', 'true']:
                                            if lastsvtype == 'TR':
                                                co = 'TR'
                                            else:
                                                co = 'Type2'
                                            f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                            manualfile.write(f)
                                        else:
                                            output(tmp_record, overjuge, 'false', real, insdelfile,
                                                   refname, numall, manualfile,lastsvtype,trfile)
                                        if real == len(numall):
                                            finaljuge = 'false'
                                            truestore.append(finaljuge)
                                        else:
                                            overjuge = 'false'
                                            truestore.append(overjuge)
                                    else:

                                        overjuge = 'false'
                                        truestore.append(overjuge)

                                    tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                  svtype, overlap, overwho, sub, \
                                                  refori_s, \
                                                  f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                  qryori_s, \
                                                  f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                    lastsvtype = 'INS'

                                else:
                                    svtype = 'DEL'
                                    overfind = 'no'
                                    overlap = '1'
                                    overwho = 'ref'
                                    sub = '0'
                                    refsv_S = str(ref1)
                                    refsv_E = str(ref2)
                                    qrysv_S = str(qry2)
                                    qrysv_E = str(qry1)
                                    f_ref_e1 = 'uncertain'
                                    f_qry_e1 = 'uncertain'
                                    f_ref_s2 = 'uncertain'
                                    f_ref_e2 = ref_oppe
                                    f_qry_s2 = 'uncertain'
                                    f_qry_e2 = qry_oppe
                                    # n1 = numall[numall.index(real)]
                                    # n2 = numall[numall.index(real) - 1]
                                    # manualfile.write(refname + '\t' + str(n1) + ',' + str(n2)+'\t'+svtype + '\n')
                                    if real >= 3:
                                        if truestore[-2:] == ['true', 'true']:
                                            if lastsvtype == 'TR':
                                                co = 'TR'
                                            else:
                                                co = 'Type2'
                                            f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                            manualfile.write(f)
                                        else:
                                            output(tmp_record, overjuge, 'true', real, insdelfile,
                                               refname, numall, manualfile, lastsvtype,trfile)
                                        if real == len(numall):
                                            finaljuge = 'true'
                                            truestore.append(finaljuge)
                                        else:
                                            overjuge = 'true'
                                            truestore.append(overjuge)
                                            
                                    else:

                                        overjuge = 'true'
                                        truestore.append(overjuge)

                                    tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                  svtype, overlap, overwho, sub, \
                                                  refori_s, \
                                                  f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                  qryori_s, \
                                                  f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                    lastsvtype = 'DEL'

                            else:#
                                #type2-DEL:
                                if qry1 > qry2:#over
                                    overlen = abs(qry2 - qry1) + 1
                                    if overlen < 45:
                                        svtype = 'DEL'
                                        overfind = 'yes'
                                        overlap = '0'
                                        overwho = '0'
                                        sub = '0'
                                        refsv_S = str(ref1 + 1)
                                        refsv_E = str(ref2 + overlen-1-1 )
                                        qrysv_S = str(qry1)
                                        qrysv_E = str(qry1)
                                        f_ref_e1 = str(ref1)
                                        f_qry_e1 = str(qry1)
                                        f_ref_s2 = str(ref2 + overlen-1 )
                                        f_ref_e2 = ref_oppe
                                        f_qry_s2 = str(qry1)
                                        f_qry_e2 = qry_oppe
                                        newref1 = str(int(refsv_E) + 1)
                                        newqry1 = qrysv_E
                                        if real >= 3:
                                            if truestore[-2:] == ['true', 'true']:
                                                if lastsvtype == 'TR':
                                                    co = 'TR'
                                                else:
                                                    co = 'Type2'
                                                f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                                manualfile.write(f)
                                            else:
                                                output(tmp_record, overjuge, 'false', real, insdelfile,
                                                   refname, numall, manualfile, lastsvtype,trfile)
                                            if real == len(numall):
                                                finaljuge = 'false'
                                                truestore.append(finaljuge)
                                            else:
                                                overjuge = 'false'
                                                truestore.append(overjuge)
                                        else:

                                            overjuge = 'false'
                                            truestore.append(overjuge)

                                        tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                      svtype, overlap, overwho, sub, \
                                                      refori_s, \
                                                      f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                      qryori_s, \
                                                      f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                        lastsvtype = 'DEL'

                                    else:
                                        svtype = 'DEL'
                                        overfind = 'no'
                                        overlap = '1'
                                        overwho = 'ref'
                                        sub = '0'
                                        refsv_S = str(ref1)
                                        refsv_E = str(ref2)
                                        qrysv_S = str(qry2)
                                        qrysv_E = str(qry1)
                                        f_ref_e1 = 'uncertain'
                                        f_qry_e1 = 'uncertain'
                                        f_ref_s2 = 'uncertain'
                                        f_ref_e2 = ref_oppe
                                        f_qry_s2 = 'uncertain'
                                        f_qry_e2 = qry_oppe
                                        if real >= 3:
                                            if truestore[-2:] == ['true', 'true']:
                                                if lastsvtype == 'TR':
                                                    co = 'TR'
                                                else:
                                                    co = 'Type2'
                                                f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                                manualfile.write(f)
                                            else:
                                                output(tmp_record, overjuge, 'true', real, insdelfile,
                                                   refname,numall, manualfile, lastsvtype,trfile)
                                            if real == len(numall):
                                                finaljuge = 'true'
                                                truestore.append(finaljuge)
                                            else:
                                                overjuge = 'true'
                                                truestore.append(overjuge)
                                        else:

                                            overjuge = 'true'
                                            truestore.append(overjuge)

                                        tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                      svtype, overlap, overwho, sub, \
                                                      refori_s, \
                                                      f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                      qryori_s, \
                                                      f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                        lastsvtype = 'DEL'

                                # type1-DEL
                                if qry1 == qry2:#type1-DEL
                                    svtype = 'DEL'
                                    overfind = 'no'
                                    overlap = '0'
                                    overwho = '0'
                                    sub = '0'
                                    refsv_S = str(ref1 + 1)
                                    refsv_E = str(ref2 - 1)
                                    qrysv_S = str(qry1)
                                    qrysv_E = str(qry1)
                                    f_ref_e1 = str(ref1)
                                    f_qry_e1 = str(qry1)
                                    f_ref_s2 = str(ref2)
                                    f_ref_e2 = ref_oppe
                                    f_qry_s2 = str(qry1)
                                    f_qry_e2 = qry_oppe
                                    if real >= 3:
                                        if truestore[-2:] == ['true', 'true']:
                                            if lastsvtype == 'TR':
                                                co = 'TR'
                                            else:
                                                co = 'Type2'
                                            f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                            manualfile.write(f)
                                        else:
                                            output(tmp_record, overjuge, 'false', real, insdelfile, refname,
                                               numall, manualfile, lastsvtype,trfile)
                                        if real == len(numall):
                                            finaljuge = 'false'
                                            truestore.append(finaljuge)
                                        else:
                                            overjuge = 'false'
                                            truestore.append(overjuge)
                                    else:

                                        overjuge = 'false'
                                        truestore.append(overjuge)
                                    tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                  svtype, overlap, overwho, sub, \
                                                  refori_s, \
                                                  f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                  qryori_s, \
                                                  f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                    lastsvtype = 'DEL'

                                ##substitution
                                if qry1 < qry2:
                                    #sub
                                    svtype = 'INS+DEL'
                                    overfind = 'no'
                                    overlap = '0'
                                    overwho = '0'
                                    sub = '1'
                                    refsv_S = str(ref1 + 1)
                                    refsv_E = str(ref2 - 1)
                                    qrysv_S = str(qry1+1)
                                    qrysv_E = str(qry2 -1)
                                    f_ref_e1 = str(ref1)
                                    f_qry_e1 = str(qry1)
                                    f_ref_s2 = str(ref2)
                                    f_ref_e2 = ref_oppe
                                    f_qry_s2 = str(qry1)
                                    f_qry_e2 = qry_oppe
                                    if real >= 3:
                                        if truestore[-2:] == ['true', 'true']:
                                            if lastsvtype == 'TR':
                                                co = 'TR'
                                            else:
                                                co = 'Type2'
                                            f = refname + '\t' + str(numall[real - 2 - 1]) +','+ str(numall[real - 1
                                                                                                             - 1]) +'\t'+ co + '\t' + '\n'
                                            manualfile.write(f)
                                        else:
                                            output(tmp_record,overjuge,'false',real,insdelfile, refname,
                                               numall,manualfile, lastsvtype,trfile)
                                        if real == len(numall):
                                            finaljuge = 'false'
                                            truestore.append(finaljuge)
                                        else:
                                            overjuge = 'false'
                                            truestore.append(overjuge)
                                    else:

                                        overjuge = 'false'
                                        truestore.append(overjuge)

                                    tmp_record = [refchr, refsv_S, refsv_E, qrychr, qrysv_S, qrysv_E, \
                                                  svtype, overlap, overwho, sub, \
                                                  refori_s, \
                                                  f_ref_s1, f_ref_e1, f_ref_s2, f_ref_e2, \
                                                  qryori_s, \
                                                  f_qry_s1, f_qry_e1, f_qry_s2, f_qry_e2]
                                    lastsvtype = 'INS+DEL'

                            if overfind == 'yes':
                                if lastsvtype == 'INS':
                                    f_qry_s1 = newqry1
                                    f_ref_s1 = newref1
                                if lastsvtype == 'DEL':
                                    f_ref_s1 = newref1
                                    f_qry_s1 = newqry1
                            else:
                                f_ref_s1 = ref_opps
                                f_qry_s1 = qry_opps
                            ref1 = int(ref_oppe)
                            qry1 = int(qry_oppe)
                            copy = 1
        # print(numall.index(real))
        try:
            if real == 2:
                # print(numall)
                output(tmp_record, overjuge, 'false', real+1, insdelfile, refname, numall, manualfile, lastsvtype,trfile)
            else:
                if truestore[-2:] == ['true', 'true']:
                    if lastsvtype == 'TR':
                        co = 'TR'
                    else:
                        co = 'Type2'
                    f = refname + '\t' +str(numall[real+1 - 2 - 1])+','+ str(numall[real+1 - 1 - 1]) +'\t'+ co + '\t'\
                        + \
                        '\n'
                    manualfile.write(f)
                else:
                    # print('1')
                    output(tmp_record, finaljuge,'false',real+1, insdelfile, refname, numall, manualfile,
                           lastsvtype,trfile)
        except NameError:
            # print('1')
            pass
        mcoordfile.close()
        manualfile.close()
        insdelfile.close()
        trfile.close()

    if model == 'batchmanual':
        manualfile = open(manualname)
        standfile = open(standname, 'w')

        tmp = []
        l = 0
        f1 = ''
        while True:
            lines = manualfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    refchr = cut[0]
                    refsvs = cut[1]
                    refsve = cut[2]
                    qrychr = cut[3]
                    qrysvs = cut[4]
                    qrysve = cut[5]
                    svtype = cut[6]
                    reforis = cut[7]
                    reffs1 = cut[8]
                    reffe1 = cut[9]
                    reffs2 = cut[10]
                    reffe2 = cut[11]
                    qryoris = cut[12]
                    qryfs1 = cut[13]
                    qryfe1 = cut[14]
                    qryfs2 = cut[15]
                    qryfe2 = cut[16]
                    l += 1
                    if reforis not in tmp:
                        tmp.append(reforis)
                        if l == 1:
                            pass
                        else:
                            f1[8] = '0'
                            f1[9] = '0'
                            f1[13] = '0'
                            f1[14] = '0'
                            standfile.write('\t'.join(f1) + '\n')
                        f1 = cut
                        n = 1
                    else:
                        n += 1
                        if n > 2:
                            f1[8] = '0'
                            f1[9] = '0'
                            f1[13] = '0'
                            f1[14] = '0'
                        f1[11] = reffe1
                        f1[16] = qryfe1
                        standfile.write('\t'.join(f1)+'\n')
                        f1 = cut
        f1[8] = '0'
        f1[9] = '0'
        f1[13] = '0'
        f1[14] = '0'
        standfile.write('\t'.join(f1) + '\n')
        standfile.close()
        manualfile.close()

    if model == 'turntovcf':
        # import numpy as np
        from pyfaidx import Fasta
        standfile = open(standname)
        vcffile = open(vcfname,'w')
        refseqfile = Fasta(refname)
        qryseqfile = Fasta(qryname)
        correct = open(correctname,'w')
        ####add header
        ##orignial header :
        header_lines = [
            "##fileformat=VCFv4.2",
            "##source=Genetic Variants Benchmark",
            "##reference=SL4.0",
            "##sample=VF36",
            "##contig=<ID=1,length=90863682>",
            "##contig=<ID=2,length=53473368>",
            "##contig=<ID=3,length=65298490>",
            "##contig=<ID=4,length=64459972>",
            "##contig=<ID=5,length=65269487>",
            "##contig=<ID=6,length=47258699>",
            "##contig=<ID=7,length=67883646>",
            "##contig=<ID=8,length=63995357>",
            "##contig=<ID=9,length=68513564>",
            "##contig=<ID=10,length=64792705>",
            "##contig=<ID=11,length=54379777>",
            "##contig=<ID=12,length=66688036>",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">",
            "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of variant\">",
            "##INFO=<ID=INS,Number=1,Type=String,Description=\"This type represents an Insertion.\">",
            "##INFO=<ID=DEL,Number=1,Type=String,Description=\"This type represents a Deletion\">",
            "##INFO=<ID=INS+DEL,Number=1,Type=String,Description=\"This type represents a Substitute\">",
            "##INFO=<ID=INV,Number=1,Type=String,Description=\"This type represents a Inversion\">",
            "##INFO=<ID=SNP,Number=1,Type=String,Description=\"This type represents a Single nucleotide "
            "polymorphism\">",
            # "##INFO=<ID=Indel_ins,Number=1,Type=String,Description=\"This type represents a Inversion\">",
            "##INFO=<ID=REFEND,Number=1,Type=Integer,Description=\"End position of the variant in the reference "
            "genome\">",
            "##INFO=<ID=REFLEN,Number=1,Type=Integer,Description=\"Length of the variant on the reference genome\">",
            "##INFO=<ID=ALTLEN,Number=1,Type=Integer,Description=\"Length of the variant on the sample genome\">",
            "##INFO=<ID=ALTPOS,Number=1,Type=String,Description=\"Range of the variant on the sample genome\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        ]
        headerori = '\n'.join(header_lines)
        vcffile.write(headerori+'\n')

        QUAL = '.'
        FILTER = 'PASS'
        SAMPLE = '1|1'
        FORMAT = 'GT'
        n = 0
        while True:
            lines = standfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    refchr = cut[0]
                    refsvs = cut[1]
                    refsve = cut[2]
                    qrychr = cut[3]
                    qrysvs = cut[4]
                    qrysve = cut[5]
                    svtype = cut[6]
                    sub = cut[7]
                    reforis = cut[8]
                    reffs1 = cut[9]
                    reffe1 = cut[10]
                    reffs2 = cut[11]
                    reffe2 = cut[12]
                    qryoris = cut[13]
                    qryfs1 = cut[14]
                    qryfe1 = cut[15]
                    qryfs2 = cut[16]
                    qryfe2 = cut[17]
                    #
                    CHROM = refchr
                    POS = refsvs

                    reflen = int(refsve) - int(refsvs) + 1
                    qrylen = int(qrysve) - int(qrysvs) + 1
                    if svtype == 'INS+DEL':
                        # if qrylen == 1:
                        #     svtype = 'DEL'
                        # if reflen == 1:
                        #     svtype = 'INS'
                        if reflen == 0:
                            svtype = 'INS'
                            refsvs = refsve
                            reffs2 = reffe1
                            sub = '0'
                        if qrylen == 0:
                            svtype = 'DEL'
                            qrysvs = qrysve
                            qryfs2 = qryfe1
                            sub = '0'
                        if qryfe1 == qryfs2:
                            qryfs2 = str(int(qrysve)+1)

                    if int(refsve) < int(refsvs) or int(qrysve) < int(qrysvs):
                        # vcffile.write(line+'\n')
                        n += 1
                    REF = refseqfile[refchr][int(refsvs) - 1:int(refsve)].seq
                    ALT = qryseqfile[qrychr][int(qrysvs) - 1:int(qrysve)].seq
                    ID = svtype+'-'+POS
                    TYPE = 'TYPE='+svtype
                    REFEND = 'REFEND='+refsve
                    # if svtype == 'INS+DEL' or svtype == 'DEL' or svtype == 'Indel-del' or svtype == 'INV':
                    #     length = int(refsve) - int(refsvs) + 1
                    # else:
                    #     length = int(qrysve) - int(qrysve) + 1
                    REFLEN = 'REFLEN='+str(int(refsve) - int(refsvs) + 1)
                    ALTLEN = 'ALTLEN='+str(int(qrysve) - int(qrysvs) + 1)
                    ALTPOS = "ALTPOS="+qrychr+':'+qrysvs+'-'+qrysve
                    INFOpre = [TYPE,REFEND,REFLEN,ALTPOS,ALTLEN]
                    INFO = ';'.join(INFOpre)
                    finalpre = [CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE]
                    vcffile.write('\t'.join(finalpre)+'\n')
                    cofinal = [refchr,refsvs,refsve,qrychr,qrysvs,qrysve,svtype,sub,reforis,reffs1,reffe1,reffs2,
                               reffe2,qryoris,qryfs1,qryfe1,qryfs2,qryfe2]
                    cofinal1 = '\t'.join(cofinal)
                    correct.write(cofinal1+'\n')
        # print(n)
        standfile.close()
        vcffile.close()
        correct.close()

    if model == 'mcoordstrue':
        mcoordfile = open(mcoordname)
        outputfile = open(outputname, 'w')
        while True:
            lines = mcoordfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    qrys = cut[0]
                    qrye = cut[1]
                    refs = cut[2]
                    refe = cut[3]
                    qrylen = cut[4]
                    reflen = cut[5]
                    similar = cut[6]
                    qryinfo = cut[11].split(':')[1].split('-')
                    qryalls = int(qryinfo[0])
                    refinfo = cut[12].split(':')[1].split('-')
                    refalls = int(refinfo[0])
                    realqrys = qryalls + int(qrys) -1
                    realqrye = qryalls + int(qrye) -1
                    realrefs = refalls + int(refs) - 1
                    realrefe = refalls + int(refe) - 1
                    final = [qrys,qrye,str(realqrys),str(realqrye),'|',refs,refe,str(realrefs),str(realrefe),'|',
                             qrylen,reflen,'|',similar,str(qryalls),
                             str(refalls)]
                    f = '\t'.join(final)
                    outputfile.write(f+'\n')
        mcoordfile.close()
        outputfile.close()
                    
    if model == 'svformatchange':
        inputfile = open(input)
        outputfile = open(outputname, 'w')

        def determine_sv_type(s):

            sv_types = ['INS', 'DEL', 'DUP', 'INV']

            matches = [sv for sv in sv_types if sv in s]

            if len(matches) == 1:
                if '01' in s:
                    svtype = 'Complex:'+matches[0]+'-DUP'
                else:
                    svtype = matches[0]
            elif len(matches) > 1:
                svtype = 'Complex:'+'-'.join(matches)
            else:
                svtype = 'DEL'
            return svtype


        def replace_callers(input_string):

            replacements = {
                'svmu': 'SVMU',
                'svision': 'SVision',
                'ngmlr': 'NGMLR',
                'anchrowave': 'AnchorWave',
                'anchorwave': 'AnchorWave',
                'syri': 'SyRI',
                'pbsv': 'PBSV',
                'pbmm2': 'pbmm2',
                'assemblytics': 'Assemblytics',
                'sniffles':'Sniffles',
                'cutesv':'cuteSV',
                'svim':'SVIM'
            }
            for old, new in replacements.items():
                input_string = input_string.replace(old, new)
            return input_string


        while True:
            lines = inputfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    s = cut[1]
                    e = cut[2]
                    svtype = cut[3]
                    callers = cut[4]
                    new_svtype = determine_sv_type(svtype)
                    if new_svtype == 'INS':
                        e = s
                    if 'Complex' in new_svtype:
                        anno = new_svtype.split(':')[1]
                        new_svtype = 'Complex'
                    else:
                        anno = 'none'
                    formatted_string = replace_callers(callers)
                    out = cut[0]+'\t'+s+'\t'+e+'\t'+new_svtype+'\t'+anno+'\t'+formatted_string+'\n'
                    outputfile.write(out)
        inputfile.close()
        outputfile.close()

    if model == 'bedtoregion':
        inputfile = open(inputname)
        outputfile = open(outputname, 'w')

        def can_merge(region1, region2):
            start1, end1 = region1
            start2, end2 = region2
            return max(start1, start2) <= min(end1, end2)
        l=0
        while True:
            lines = inputfile.readlines(10000)
            if not lines:
                break
            else:
                for line in lines:
                    line = line.rstrip()
                    cut = line.strip().split('\t')
                    chr = cut[0]
                    s = int(cut[1])
                    e = int(cut[2])
                    l+=1
                    if chr == 'Chr01' and cut[4] == '60':
                        if l == 1:
                            region1 =(s,e)
                        else:
                            now_region = (s,e)
                            merge = can_merge(region1, now_region)
                            

