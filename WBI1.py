import numpy as np
import operator
def wczytaj (plik):
    p = open(plik, "r+")
    p.readline()
    s=""
    #p.truncate(0)
    for line in p:
        s=s+line
    licz=0

    s1=list(s)
    print ("aaa", s1[licz])
    while (licz<len(s1)):
        if (s1[licz]=="A"):
            s1[licz]="U"
        elif (s1[licz]=="G"):
            s1[licz]="C"
        elif (s1[licz]=="C"):
            s1[licz]="G"
        elif (s1[licz] == "T"):
            s1[licz] = "A"
        licz=licz+1;

    p1=open("kolagen.fasta", "w")
    p1.write("".join(s1))
    p.close()
    p1.close()
    print("abc")

def traf (sym1, sym2):
    if sym1 == sym2:
        return 3
    else:
        return -3
def needleman_wunsch (s1,s2):
    kara = -1

    macierz = np.zeros((len(s1)+1,len(s2)+1),dtype='i,i,i')
    seq1=list(s1)
    seq2=list(s2)
    for i in range(len(s2)+1):
        macierz[0][i]=(i*-1,0,i)
    for j in range(len(s1)+1):
        macierz[j][0]=(j*-1,j,0)
    for ii in range (1,len(s1)+1):
        for jj in range (1, len(s2)+1):

            lis=[(traf(seq1[ii - 1], seq2[jj - 1]) + macierz[ii - 1][jj - 1][0], ii - 1, jj - 1), (macierz[ii - 1][jj][0] + kara, ii - 1, jj), (macierz[ii][jj - 1][0] + kara, ii, jj - 1)]
            macierz[ii][jj]=max(lis, key=operator.itemgetter(0))

    print(macierz)
    print(seq1)
    print(seq2)
    start=macierz[len(s1)][len(s2)]

    strt=""
    strt2=s1[start[1]]
    print("start", start[1])
    flag=0

    while (start[1]!=0 and start[2]!=0):
        elem=start
        start=macierz[start[1]][start[2]]
        if start[1]!=elem[1] and start[2]==elem[2] :
            strt+="_"
            strt+=s2[elem[1]]
            strt2 += s1[start[1]]
        elif start[1]!=elem[1] and start[2]!=elem[2] :
            strt+=s2[elem[1]]
            strt2+=s1[start[1]]
        elif start[1]==elem[1] and start[2]!=elem[2] :
            strt2 += "_"

    strt += s2[start[1]]
    strt=strt[::-1]
    strt2=strt2[::-1]
    print ("uuu", strt)
    print("uuu2", strt2)
def smith_waterman(s1,s2):
    kara=-2
    macierz = np.zeros((len(s1) + 1, len(s2) + 1))
    seq1 = list(s1)
    seq2 = list(s2)
    for ii in range (1,len(s1)+1):
        for jj in range (1, len(s2)+1):

            lis=[0, (traf(seq1[ii - 1], seq2[jj - 1]) + macierz[ii - 1][jj - 1]), (macierz[ii - 1][jj] + kara), (macierz[ii][jj - 1] + kara)]

            macierz[ii][jj]=max(lis)
    macierz=macierz.transpose()
    max_value=macierz.max()
    #macierz[0][0]=13
    print(macierz)
    coordinates=np.where(macierz == max_value)
    pkt_startu = np.zeros((len(coordinates[0])), dtype='i,i')
    str1=[""] * len(coordinates[0])
    str2 = [""] * len(coordinates[0])
    for i in range (len(coordinates[0])):
        pkt_startu[i]=(coordinates[0][i],coordinates[1][i])
        licz1=0
        licz2=0
        eee = pkt_startu[i][1]
        while (macierz[pkt_startu[i][0]-licz1][pkt_startu[i][1]-licz2]>0 and (pkt_startu[i][0]-licz1>=0 and pkt_startu[i][0]-licz2)>=0 ): #czy wieksze czy rowne?

            temp=macierz[pkt_startu[i][0]-licz2][pkt_startu[i][1]-licz1]
            #zamiast licz+= przek= a[-1][-1], lew=a[0][-1] gora [-1][0] if przek>gora i przek > lewo to licz1, licz2-=1 elif lew>przek i lew>gora to licz1-=1 licz2=licz2, licz1=licz1, licz2-=1
            przekatna = macierz[pkt_startu[i][0]-licz2-1][pkt_startu[i][1]-licz1-1]
            lewo=macierz[pkt_startu[i][0]-licz2][pkt_startu[i][1]-licz1-1]
            gora = macierz[pkt_startu[i][0]-licz2-1][pkt_startu[i][1]-licz1]
            if (przekatna >= lewo and przekatna >=gora) or (gora <=temp and lewo <=temp):
                str1[i] = str1[i] + s1[pkt_startu[i][1] - licz1 - 1]  # jeszcze podłogi
                str2[i] = str2[i] + s2[pkt_startu[i][1] - licz2]
                licz1+=1
                licz2+=1
            elif (przekatna <= lewo and lewo >=gora):
                str2[i] = str2[i] + "_"
                str1[i] = str1[i] + s2[pkt_startu[i][1] - licz1]
                licz1+=1
                #str1[i] = str1[i] + "_"
            elif (gora > przekatna and gora >= lewo and gora > temp):
                str1[i] = str1[i] +"_"
                str2[i]=str2[i] + s2[pkt_startu[i][1] - licz2]
                licz2 += 1
    str1=''.join(str1)
    str1 = str1[::-1]
    str2 = ''.join(str2)
    str2 = str2[::-1]
    print(pkt_startu)
    print(str1)
    print(str2)

wczytaj("kol.fasta")
#needleman_wunsch("GATTACA","GCATGCU")
smith_waterman("TGTTACGG","GGTTGACTA")
