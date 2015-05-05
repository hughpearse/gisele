#!/usr/bin/python 
from math import *


def needlemanwunsch(strA, strB):
    strG = ''.join(strA)
    strR = ''.join(strB)
    matrix = []
    path = []
    row = len(strR) 
    col = len(strG)
    strR = "^"+strR
    strG = "^"+strG
    for i in range(row+1):
        matrix.append([0]*(col+1))
        path.append(["N"]*(col+1))
    
    def print_matrix(matrix):
        print '\t'+('\t'.join(map(str,list(strG))))
        i = 0
        for line in matrix:
            print strR[i]+"\t"+('\t'.join(map(str,line)))
            i +=1
    
    # print_matrix(matrix)
    indelValue = -1
    matchValue = 2
    for i in range(1,row+1):
        for j in range(1,col+1):
            # penalty map
            from_left = matrix[i][j-1] + indelValue
            from_top = matrix[i-1][j] + indelValue
            if strR[i]==strG[j]:
                from_diag = matrix[i-1][j-1] + matchValue
            else:
                from_diag = matrix[i-1][j-1] + indelValue
    
            matrix[i][j]= max(from_left,from_top,from_diag)
            # path map
            if matrix[i][j]==from_left:
                path[i][j]="-"
            elif matrix[i][j]==from_top:
                path[i][j] = "|"
            elif matrix[i][j] == from_diag:
                path[i][j] = "M"
            else:
                pass
    
            if matrix[i][j]<0:
                matrix[i][j]=0
    pass
    #print_matrix(matrix)
    #print
    #print_matrix(path)
    
    iRow = len(matrix)-1
    jCol = len(matrix[0])-1
    
    while iRow>=0:
        maxPnltyVaue = max(matrix[iRow])
        while jCol>=0:
            if matrix[iRow][jCol]==maxPnltyVaue:
                ipair = iRow
                jpair = jCol
                reportR=[]
                reportG=[]
                while (1):
                    if ipair ==0 and jpair==0:
                        break
                    # else:
                    if path[ipair][jpair]=="M":
                        reportR.append(strR[ipair])
                        reportG.append(strG[jpair])
                        ipair -= 1
                        jpair -= 1
                    elif path[ipair][jpair]=="-":
                        # reportR.append(strR[ipair])
                        reportR.append("-")
                        reportG.append(strG[jpair])
                        # ipair -= 1
                        jpair -= 1
                    elif path[ipair][jpair]=="|":
                        reportR.append(strR[ipair])
                        reportG.append("-")
                        ipair -= 1
                    elif path[ipair][jpair]=="N":
                        if ipair > 0:
                            reportR.append(strR[ipair])
                            ipair -=1
                        if jpair > 0:
                            reportG.append(strG[jpair])
                            jpair -= 1
                s1 = "".join(reportR)
                s2 = "".join(reportG)
                
            jCol -=1
        iRow -=1
    
    s3 = ""
    i=0
    while i < max(len(s1), len(s2)):
        try:
            if(s1[i] == s2[i]):
                s3 += s1[i]
            else:
                s3 += "-"
        except IndexError:
            s3 += "-"
        i+=1
    
    return list("".join(reversed(s3)))

#strG = "isaac went to the shop"
strG = ["e", "\xad", "\xbf", "\xef"]
#strR = "ezra went shopping on his own"
#strR = "e\xad\xbe\xef"
strR = ["e", "\xad", "\xdf", "\xef"]
gms = needlemanwunsch(strG, strR)
print gms