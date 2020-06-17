import pysam 
import os 
import re
import numpy as np
import pandas as pd 
import time
from pysam import VariantFile
from sklearn.cluster import KMeans
import sklearn.cluster
pd.options.mode.chained_assignment = None


def CallVCF(CandidateDf, VCFFile, refFile, bamFile, CellName):
    # transfer candidate format to vcf format 
    ## Select Proper SNVs 
    # Exact 0/0 position will not appeared in the vcf
    Condition1 = CandidateDf['Genotype'].isin(['0/1', '1/1'])
    Condition2 = CandidateDf['Decision'].isin(['PCRError', 'PCRErrorLow'])
    Condition2_1 = CandidateDf['Decision'].isin(['PCRError'])
    Condition2_2 = CandidateDf['Decision'].isin(['PCRErrorLow'])
    Condition3 = CandidateDf['AlleleDropOut'].isin([1])
    Condition4 = CandidateDf['Vague'].isin(['Yes'])
    Condition5 = CandidateDf['DataRange'].isin(['HomoHeteroSNV'])
    Condition6 = CandidateDf['PassCode'].isin(['ClusterSNV'])
    ## make a VCF head 
    # Info items include: NS, DP, AF, AA
    # Filter Items include: AB (Amplicon Bias), PCRError (PCR Error), ADO (Allele drop out)
    # Format Items include: GT (0/1, 1/1, 0/0) , L0 (Homozygous Likelihood), L1 (Amplicon Bias Likelihood), L2 (Heterozygous Likelihood), L3 (PCR error Likelihood), DP, AF
    with open(VCFFile, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n##fileDate=%s\n##source=NBJUDGE\n")
        vcf.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n")
        vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        vcf.write("##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n")
        vcf.write("##FILTER=<ID=AB,Description=\"Amplicon Bias\">\n")
        vcf.write("##FILTER=<ID=PCRError,Description=\"Suspected PCR Error\">\n")
        vcf.write("##FILTER=<ID=ADO,Description=\"Suspected Allele Drop-out Error\">\n")
        vcf.write("##FILTER=<ID=NE,Description=\"Not enough homozygous or heterozygous SNVs in the neighborhood\">\n")
        vcf.write("##FILTER=<ID=CC,Description=\"SNV located in SNV cluster\">\n")
        vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        vcf.write("##FORMAT=<ID=L0,Number=1,Type=Float,Description=\"Homozygous likelihood\">\n")
        vcf.write("##FORMAT=<ID=L1,Number=1,Type=Float,Description=\"Amplicon Bias Likelihood\">\n")
        vcf.write("##FORMAT=<ID=L2,Number=1,Type=Float,Description=\"Heterozygous likelihood\">\n")
        vcf.write("##FORMAT=<ID=L3,Number=1,Type=Float,Description=\"PCR Error likelihood\">\n")
        vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        vcf.write("##FORMAT=<ID=AF,Number=4,Type=Integer,Description=\"Top 4 Base count\"\n")
        for chrom in set(CandidateDf['Chrom']):
            vcf.write('##contig=<ID=%s>\n' % chrom)
        vcf.write("##reference=file://%s\n" % refFile)
        vcf.write("##bamFile=file://%s\n" % bamFile)
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % CellName)
        vcf.flush()
        # AddLabel For each SNV in the CandidateDf_VCF 
        CandidateDf['AF'] = CandidateDf['AltCount'] / CandidateDf['DP'] 
        CandidateDf['INFO'] = 'NS=1;DP=' + CandidateDf['DP'].astype(str) + ";AF=" + CandidateDf['AF'].astype(str) + ";AA=" + CandidateDf['Alt']
        CandidateDf['FILTER'] = ""
        CandidateDf.loc[Condition2_1, 'FILTER'] += 'AB'
        CandidateDf.loc[Condition2_2, 'FILTER'] += 'PCRError'
        Condition_TMP = CandidateDf['FILTER'].isin([""])
        CandidateDf.loc[Condition3&Condition4&Condition_TMP, 'FILTER'] += 'ADO'
        CandidateDf.loc[Condition3&Condition4&(~Condition_TMP), 'FILTER'] += ',ADO'
        Condition_TMP = CandidateDf['FILTER'].isin([""])
        CandidateDf.loc[(~Condition5)&Condition_TMP, 'FILTER'] += 'NE'
        CandidateDf.loc[(~Condition5)&(~Condition_TMP), 'FILTER'] += ',NE'
        Condition_TMP = CandidateDf['FILTER'].isin([""])
        CandidateDf.loc[(Condition6)&Condition_TMP, 'FILTER'] += 'CC'
        CandidateDf.loc[(Condition6)&(~Condition_TMP), 'FILTER'] += ',CC'
        Condition_TMP = CandidateDf['FILTER'].isin([""])
        CandidateDf.loc[Condition_TMP, 'FILTER'] = 'PASS'
        CandidateDf['FORMAT_Label'] = 'GT:L0:L1:L2:L3:DP:AF'
        CandidateDf['FORMAT'] = CandidateDf['Genotype'] + ":" + CandidateDf['HomoZygous'].astype(str) + ":" + CandidateDf['HeteroZygous'].astype(str) + ":" + CandidateDf['PCRError'].astype(str) + ":" + CandidateDf['PCRErrorLow'].astype(str) + ":" + CandidateDf['DP'].astype(str) + ":" + CandidateDf['Base1'].astype(str) + "," + CandidateDf['Base2'].astype(str) + "," + CandidateDf['Base3'].astype(str) + "," + CandidateDf['Base4'].astype(str)
        CandidateDf_VCF = CandidateDf[Condition1 | Condition2 | Condition3 | Condition4]
        CandidateDf_VCF['VCFInfo'] = CandidateDf_VCF['Chrom'] + "\t" + CandidateDf_VCF['Pos'].astype(str) + "\t.\t" + CandidateDf_VCF['Ref'] + "\t" + CandidateDf_VCF['Alt'] + "\t.\t" + CandidateDf_VCF['FILTER'] + '\t' + CandidateDf_VCF['INFO'] + "\t" + CandidateDf_VCF['FORMAT_Label'] + "\t" + CandidateDf_VCF['FORMAT'] 
        
        for I in CandidateDf_VCF.index:
            Record = CandidateDf_VCF.loc[I, 'VCFInfo']
            vcf.write(Record + "\n")
            vcf.flush()
def WeightKernel(Distance, mu=np.log(3)/10):
    # return weight vector 
    fx = 2 * np.exp(-mu * Distance) / (1 + np.exp(-mu * Distance))
    return(fx)
    

def CalculateP_A(SubDf, Position, mu=np.log(3)/10):
    # set kernal 
    SubDf['Distance'] = np.abs(np.array(SubDf['Pos'])-Position) / 1000.0
    WorkDf = SubDf.loc[SubDf['Distance']!=0]
    Weight = WeightKernel(np.array(WorkDf['Distance'], dtype=float), mu)
    ProbDict = {}
    for Base_Index in range(1,5):
        BaseName = 'Base%s' % Base_Index
        RateName = 'BRate%s' % Base_Index
        WorkDf[RateName] = WorkDf[BaseName] / WorkDf['DP']
        Ex = np.sum(np.array(WorkDf[RateName], dtype=float) * Weight) / np.sum(Weight)
        Varx = np.sqrt(np.sum(Weight * (np.array(WorkDf[RateName], dtype=float) - Ex)**2) / np.sum(Weight))
        ProbDict[BaseName] = [Ex, Varx]
    return(ProbDict)
    
def DfTailor(SubDf1, SubDf2, I):
    # gain Sub dataframe with tha same dim (decided by the smaller Dataframe)
    # Tail the larger dataframe around the candidate Point
    # I is the position value of candidate SNV point
    # discard the current position
    SubDf1 = SubDf1.loc[SubDf1['Pos']!=I]
    SubDf2 = SubDf2.loc[SubDf2['Pos']!=I]
    L1 = SubDf1.shape[0]
    L2 = SubDf2.shape[0]
    SubDf1['Distance'] = np.abs(np.array(SubDf1['Pos'], dtype=int) - I)
    SubDf2['Distance'] = np.abs(np.array(SubDf2['Pos'], dtype=int) - I)
    S_DF = SubDf1
    L_DF = SubDf2
    if L2 < L1:
        S_DF = SubDf2
        L_DF = SubDf1
    # large Df should have Position that is similar to small Df
    LocList = []
    S_DF.sort_values("Distance", inplace=True)
    for Site in S_DF.index:
        tmpPos = S_DF.loc[Site, 'Pos']
        L_DF['tmpDistance'] = np.abs(np.array(L_DF['Pos'], dtype=int)-tmpPos)
        L_DF.sort_values("tmpDistance", inplace=True)
        LocAdd = 0
        for i in range(0,L_DF.shape[0]):
            if L_DF.index[i] not in LocList:
                LocList.append(L_DF.index[i])
                LocAdd += 1
                break
        if LocAdd == 0:
            S_DF = S_DF.drop(index=[Site])
    L_DF = L_DF.loc[LocList]
    ResultDF = pd.concat([S_DF, L_DF], sort=True)
    return(ResultDF)

def GetNeighbor(SubDf2, I):
    SubDf2['Distance'] = [abs(x) for x in list(SubDf2['Pos']-I)]
    SubDf2 = SubDf2.loc[SubDf2['Distance']!=0]
    SubDf2.sort_values("Distance", inplace=True)
    ResultDf = SubDf2.loc[SubDf2.index[0:10]]
    return(ResultDf)

def GetNeighbor_Normal(CountDf, I, WindowSize):
    '''
    Input CountDf & Pos data
    Return SNVs within 30000bp
    '''
    start = I - WindowSize
    end = I + WindowSize
    return(CountDf.loc[(CountDf['Pos']>=start)&(CountDf['Pos']<=end)])

def highlight(data, Index,color='darkorange'):
    '''
    Highligh the specific row
    '''
    ColorDf = pd.DataFrame('', columns=data.columns, index=data.index)
    ColorDf.loc[Index] = 'background-color: {}'.format(color)
    return(ColorDf)

def MakeShowData(InputDf, PCRDf, InterestedColumn1):
    HomoDf = InputDf.loc[InputDf['RawGT']=='Homozygous']
    HeterDf = InputDf.loc[InputDf['RawGT']=='Heterozygous']
    PCR_eg = {}
    for I in PCRDf.index:
        currentSite = int(PCRDf.loc[I, 'Pos'])
        if PCRDf.loc[I, 'DataRange'] == 'HomoHeter':
            WorkDf = GetNeighbor_Normal(InputDf, currentSite, WindowSize)
        else:
            WorkDf = pd.concat([GetNeighbor(HomoDf, currentSite), GetNeighbor(HeterDf, currentSite)], sort=True)
            WorkDf.sort_values("Pos", inplace=True)
        df = WorkDf[InterestedColumn1]
        PCR_eg[I] = df.style.apply(highlight, Index=I, axis=None)
    return(PCR_eg)

## Finish Common Functions

class EstModel():
    
    def __init__(self, bamFile, ref, CellName=None, WindowSize=30000, mu=np.log(3)/10):
        '''
        Load requriement File Names and File Path
        '''
        try:
            M = pysam.FastaFile(ref)
        except:
            print("Please identify proper reference file path")
            sys.exit(1)
        try:
            L = pysam.AlignmentFile(bamFile)
        except:
            print("Please identify proper bam file path")
            sys.exit(1)
        self.bam = pysam.AlignmentFile(bamFile)
        self.ref = pysam.FastaFile(ref)
        self.WindowSize = WindowSize
        self.mu = mu
        self.Candidate = pd.DataFrame()
        self.Shadow = None
        self.Result = None
        self.MergeResult = None
        self.ABdf = pd.DataFrame()
        self.GapList = []     # record genome region which have depth lower than 10X
        self.JUDGEMENT = None
        self.Name = os.path.basename(bamFile).split(".bam")[0]
        if not CellName:
            self.Name = CellName
    
    def MakeCandidateDf(self, chrom, Start, End):
        '''
        Convert sorted bam into candidate dataframe
        input:
            refFile: pysam.FastarFile object for reference genome
            samfile: pysam.AlignmentFile object for sorted bam file
            chrom: aim chromosome required for dataframe extraction
            Start: start position for dataframe 
            End: end position for dataframe
        '''
        self.Start = Start
        self.End = End
        if self.Start >= self.End:
            print('Too Small data range, exit!')
            sys.exit(1)
        refFile = self.ref
        samfile = self.bam
        WindowSize=self.WindowSize
        CountDf = pd.DataFrame(columns=['Chrom', 'Pos', 'Ref', 'RefCount', 'Alt', 'AltCount', 'Base1', 'Base2', 'Base3', 'Base4', 'RawRate'])
        str_ref = refFile.fetch(chrom, Start - 1, End + 1)
        my_arg = {"fastafile": refFile, "stepper": "samtools", "adjust_capq_threshold": 50, "contig": chrom,
                  "start": Start, "stop": End, "min_mapping_quality": 40, "truncate":True, "min_base_quality":17}
        CandidateDf = pd.DataFrame(columns=['Chrom', 'Pos', 'Ref', 'RefCount', 'Alt', 'AltCount', 'DP', 'Base1', 'Base2', 'Base3', 'Base4', 'RawRate'])
        partner = 0                                   # means no need for 0/0 partner in the Raw genotype decision
        partnerDf = pd.DataFrame(columns=['Chrom', 'Pos', 'DP', 'Base1', 'Base2', 'Base3', 'Base4', 'RawRate'])
        for pileup_column in samfile.pileup(**my_arg):
            pos = pileup_column.reference_pos + 1   # 1-based location
            if pos > my_arg['stop'] + self.WindowSize:
                break
            if pos < my_arg['start']-self.WindowSize:
                continue
            read_bases_list = pileup_column.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True)
            DP = 0
            RefCount = 0
            AltCount = 0
            BasePool = {'A':0,
                        'T':0,
                        'C':0,
                        'G':0
                       }
            if len(read_bases_list) >=10:
                refSeq = str_ref[pos-my_arg['start']].upper()
                if refSeq not in ['A', 'T', 'C', 'G']:
                    continue
                for B in read_bases_list:
                    base_str = B.upper()
                    if base_str in BasePool.keys():
                        BasePool[base_str] += 1
                    else:
                        BasePool[base_str] = 1
                RefCount = BasePool[refSeq]
                NonRefList = sorted([BasePool[K] for K in BasePool.keys() if K!=refSeq], reverse=True)
                AltCount = NonRefList[0]
                if AltCount >= 2: # less than max allele might come from noise
                    AltBase = sorted([(K, BasePool[K]) for K in BasePool.keys() if K!=refSeq], reverse=True, key=lambda item:item[1])[0][0]
                    DP = float(sum(NonRefList[0:3]) + RefCount)
                    BaseCountList = sorted(NonRefList[0:3] + [RefCount], reverse=True)
                    Info = [my_arg['contig'], pos, refSeq, float(RefCount), AltBase, float(AltCount), DP] + BaseCountList + [BaseCountList[1] / DP]
                    CandidateDf.loc['%s:%s' % (my_arg['contig'], pos)] = Info
                    partner += 2    # Now we need a 0/0 genotype partner 
                elif  AltCount > 0:  # proper point with Alt Allele
                    if partner >=1: # Now we need a 0/0 genotype partner 
                        DP = float(sum(NonRefList[0:3]) + RefCount)
                        BaseCountList = sorted(NonRefList[0:3] + [RefCount], reverse=True)
                        Info = [my_arg['contig'], pos, DP] + BaseCountList + [BaseCountList[1] / DP]
                        partnerDf.loc['%s:%s' % (my_arg['contig'], pos)] = Info
                        partner += -1
                    else:
                        continue
                else:
                    continue
            else:
                continue
        # return the total SNV point
        print('%s: %s candidate SNVs have been loaded !' % (time.ctime(), CandidateDf.shape[0]))
        print('%s: %s shadow base points have been loaded !' % (time.ctime(), partnerDf.shape[0]))
        self.Candidate = CandidateDf
        self.Shadow = partnerDf
        # Identify homozygous and heterozygous Roughly
        self.GetCutoffSimple()

    def GetCutoffSimple(self):
        # decrease the complex of method GetCutoff()
        try:
            self.Candidate.shape
        except:
            print("No CandidateDf exist. Run MakeCandidateDf method first")
        self.Candidate['RawGT'] = None
        CutOFFDf = self.Candidate[['Pos', 'RawRate']]
        CutOFFShadow = self.Shadow[['Pos', 'RawRate']]
        CutOFFDf.index = self.Candidate.index
        CutOFFShadow.index = self.Shadow.index
        CutOFFDf['Label'] = 'UnSet'
        MixDf = pd.concat([CutOFFDf, CutOFFShadow], axis=0, sort=True).sort_values(by=['RawRate'])
        # Make cutoff every 1 MB 
        # PosMax = list(MixDf['Pos'])[-1]
        Orig = list(CutOFFDf['Pos'])[0]-1
        estimator = sklearn.cluster.AgglomerativeClustering(2)
        # get total center Mean
        estimator.fit(pd.DataFrame(MixDf['RawRate']))
        MixDf['Label'] = estimator.labels_
        MixDf.loc[MixDf['Label']==estimator.labels_[0], 'Label'] = 'Homozygous'
        MixDf.loc[MixDf['Label']==estimator.labels_[-1], 'Label'] = 'Heterozygous'
        CutOFFDf.loc[CutOFFDf.index, ['Label']] = MixDf.loc[CutOFFDf.index, ['Label']]
        # 
        while CutOFFDf.loc[(CutOFFDf['Pos']>Orig)&(CutOFFDf['Pos']<=Orig+self.WindowSize)].shape[0] >=3:
            SubDf = CutOFFDf.loc[(CutOFFDf['Pos']>Orig)&(CutOFFDf['Pos']<=Orig+self.WindowSize)]
            SubShadow = CutOFFShadow.loc[(CutOFFShadow['Pos']>Orig)&(CutOFFShadow['Pos']<=Orig+self.WindowSize)]
            SubMix = pd.concat([SubDf, SubShadow], axis=0, sort=True).sort_values(by=['RawRate'])
            estimator.fit(pd.DataFrame(SubMix['RawRate']))
            SubMix['Label'] = estimator.labels_
            SubMix.loc[SubMix['Label']==estimator.labels_[0], 'Label'] = 'Homozygous'
            SubMix.loc[SubMix['Label']==estimator.labels_[-1], 'Label'] = 'Heterozygous'
            CutOFFDf.loc[SubDf.index, ['Label']] = SubMix.loc[SubDf.index, ['Label']]
            Orig = list(SubDf['Pos'])[-1]
        # check Null cutoffs
        self.Candidate['RawGT'] = CutOFFDf['Label']
    
    def AnnotateCountDf1(self):
        # The first annotation of candidate SNV Dataframe
        self.ABdf = pd.DataFrame(columns=['HeteroBase1', 'HeteroSE', 'HomoBase1', 'HomoSE', 'PCRBase1', 'PCRSE', 'PCRLowBase1'], index=self.Candidate.index)
        CountDf = self.Candidate
        WindowSize = self.WindowSize
        mu = self.mu
        CountDf['HomoZygous'] = 1 
        CountDf['HomoZygous_P'] = 0
        CountDf['HomoZygous_P_std'] = 0
        CountDf['HeteroZygous'] = 1
        CountDf['HeteroZygous_P'] = 0
        CountDf['HeteroZygous_P_std'] = 0
        CountDf['PCRError'] = 1
        CountDf['PCRErrorLow'] = 1
        CountDf['DataRange'] = 'AllData'
        CountDf['Decision'] = 'Unknown'
        CountDf['DecisionGT'] = 'Unknown'
        TotalHomoDf = CountDf.loc[CountDf['RawGT']=='Homozygous']
        TotalHeterDf = pd.concat([CountDf.loc[CountDf['RawGT']=='Heterozygous'], self.Shadow], axis=0, sort=True)
        for I in CountDf.index:
            currentPos = CountDf.loc[I, 'Pos']
            start = currentPos - WindowSize
            end = currentPos + WindowSize
            SubDf = CountDf.loc[(CountDf['Pos']>=start)&(CountDf['Pos']<=end)&(CountDf['Pos']!=currentPos)]
            SubShadowDf = self.Shadow.loc[(self.Shadow['Pos']>=start)&(self.Shadow['Pos']<=end)&(self.Shadow['Pos']!=currentPos)]
            SubDf_Homo = pd.concat([SubDf.loc[SubDf['RawGT']=='Homozygous'], SubShadowDf],axis=0, sort=True)
            SubDf_Heter = SubDf.loc[SubDf['RawGT']=='Heterozygous']
            DataRange = ''
            if (SubDf_Homo.shape[0] < 3) and (SubDf_Heter.shape[0] < 3):
                DataRange = 'LonelySNV'
                HomoProb = CalculateP_A(TotalHomoDf, currentPos, mu)
                HeteroProb = CalculateP_A(TotalHeterDf, currentPos, mu)
                PCRProb = CalculateP_A(DfTailor(TotalHomoDf, TotalHeterDf, currentPos), currentPos, mu)
            elif (SubDf_Homo.shape[0] < 3) and (SubDf_Heter.shape[0]>=3):
                DataRange = 'HeteroEnrichSNV'
                HomoProb = CalculateP_A(TotalHomoDf, currentPos, mu)
                HeteroProb = CalculateP_A(SubDf_Heter, currentPos, mu)
                PCRProb = CalculateP_A(DfTailor(TotalHomoDf, SubDf_Heter, currentPos), currentPos, mu)
            elif (SubDf_Homo.shape[0] >= 3) and (SubDf_Heter.shape[0]<3):
                DataRange = 'HomoEnrichSNV'
                HomoProb = CalculateP_A(SubDf_Homo, currentPos, mu)
                HeteroProb = CalculateP_A(TotalHeterDf, currentPos, mu)
                PCRProb = CalculateP_A(DfTailor(SubDf_Homo, TotalHeterDf, currentPos), currentPos, mu)
            else:
                DataRange = 'HomoHeteroSNV'
                HomoProb = CalculateP_A(SubDf_Homo, currentPos, mu)
                HeteroProb = CalculateP_A(SubDf_Heter, currentPos, mu)
                PCRProb = CalculateP_A(DfTailor(SubDf_Homo, SubDf_Heter, currentPos), currentPos, mu)
            CountDf.loc[I, 'DataRange'] = DataRange
            # fill in the std value
            CountDf.loc[I, 'HomoZygous_P'] = HomoProb['Base2'][0]
            CountDf.loc[I, 'HomoZygous_P_std'] = HomoProb['Base2'][1]
            CountDf.loc[I, 'HeteroZygous_P'] = HeteroProb['Base2'][0]
            CountDf.loc[I, 'HeteroZygous_P_std'] = HeteroProb['Base2'][1]
            self.ABdf.loc[I] = [HeteroProb['Base1'][0], HeteroProb['Base1'][1], HomoProb['Base1'][0], HomoProb['Base1'][1], PCRProb['Base1'][0], PCRProb['Base1'][1], HeteroProb['Base1'][0]-3*HeteroProb['Base1'][1]]
            for Base in ['Base1', 'Base2', 'Base3', 'Base4']:
                CountDf.loc[I, 'HomoZygous'] = CountDf.loc[I, 'HomoZygous'] + (np.log10(HomoProb[Base][0]+0.0000001)*CountDf.loc[I, Base])
                CountDf.loc[I, 'HeteroZygous'] = CountDf.loc[I, 'HeteroZygous'] + (np.log10(HeteroProb[Base][0] + 0.0000001)*CountDf.loc[I, Base])
                CountDf.loc[I, 'PCRError'] = CountDf.loc[I, 'PCRError'] + (np.log10(PCRProb[Base][0]+0.0000001)*CountDf.loc[I, Base])
            PCRLow_P1 = np.max([0.0000001, (HeteroProb["Base1"][0] + 0.0000001 - 3 * HeteroProb["Base1"][1])])
            CountDf.loc[I, 'PCRErrorLow'] = (np.log10(PCRLow_P1+0.0000001)*CountDf.loc[I, "Base1"]) + (np.log10(1-PCRLow_P1+0.0000001)*CountDf.loc[I, "Base2"])
            # Normalize the percentage data
            ProbUnit = np.array(CountDf.loc[I, ['HomoZygous', 'HeteroZygous', 'PCRError', 'PCRErrorLow']])
#             ProbUnit = ProbUnit / np.sum(ProbUnit)
            CountDf.loc[I, ['HomoZygous', 'HeteroZygous', 'PCRError', 'PCRErrorLow']] = list(ProbUnit)
            DecisionList = ['Homozygous', 'Heterozygous', 'PCRError', 'PCRErrorLow']
            Decision = DecisionList[np.where(ProbUnit==np.max(ProbUnit))[0][0]]
            CountDf.loc[I, 'Decision'] = Decision
            if Decision == 'PCRError':
                # we should make a decision which GT to choose
                if CountDf.loc[I, 'HomoZygous'] < CountDf.loc[I, 'HeteroZygous']:            # HeteroZygous have larger base2 rate than HomoZygous point
                    CountDf.loc[I, 'DecisionGT'] = 'Heterozygous'     # candidate position are more likely to HomoZygous
                else:
                    CountDf.loc[I, 'DecisionGT'] = 'Homozygous'    # else candidate position are more likely to HeteroZygous
            else:
                CountDf.loc[I, 'DecisionGT'] = CountDf.loc[I, 'Decision'] 
            if Decision == 'PCRErrorLow':
                CountDf.loc[I, 'DecisionGT'] = 'Heterozygous'
        self.Candidate = CountDf
        self.ABdf['Pos'] = CountDf['Pos']
        self.ABdf.index = CountDf.index
        self.ABdf['RawGT'] = CountDf['RawGT']
        self.ABdf['RawRate'] = CountDf['RawRate']
        self.ABdf['Decision'] = CountDf['Decision']
        self.ABdf['DecisionGT'] = CountDf['DecisionGT']
        
        
    def AnnotateCountDf2(self):
        InputDf2 = self.Candidate
        WindowSize = self.WindowSize
        # Input dataframe after the first annotation
        InputDf2['AlleleDropOut'] = 0
        for I in InputDf2.index:
            currentPos = InputDf2.loc[I, 'Pos']
            start = currentPos - WindowSize
            end = currentPos + WindowSize
            SubDf = InputDf2.loc[(InputDf2['Pos']>=start)&(InputDf2['Pos']<=end)]
            SubDf_Heter = SubDf.loc[SubDf['DecisionGT']=='Heterozygous']
            if SubDf_Heter.shape[0] == 0:    # No hetero SNV around an ADO
                InputDf2.loc[I, 'AlleleDropOut'] = 1
        self.Candidate = InputDf2


    def EstimateError(self):
        '''
        Estimate 4 types of error
        Homozygous: 
        Heterozygous:
        PCR error:
        Allele drop out:
        '''
        print('%s: Start to Calculate!' % time.ctime())
        # The first Annotation
        CountDf = self.Candidate
        WindowSize = self.WindowSize
        self.AnnotateCountDf1()
        print('%s: The first Annotation finished!' % time.ctime())
        # Now add ADO score for each point
        self.AnnotateCountDf2()
        print('%s: The second Annotation finished!' % time.ctime())
        self.Result = self.Candidate
    
    def FillError(self):
        self.Candidate['Vague'] = 'No'
        for I in self.Candidate.index:
            currentPos = self.Candidate.loc[I, 'Pos']
            SubDf = self.Candidate.loc[(self.Candidate['Pos']>currentPos-self.WindowSize/3.0)&(self.Candidate['Pos']<currentPos+self.WindowSize/3.0)]
            if (SubDf.shape[0] > 10) and (100.0 * SubDf.loc[SubDf['Decision']=='PCRError'].shape[0] / SubDf.shape[0] >= 40):
                self.Candidate.loc[SubDf.index, 'Vague'] = 'Yes'

    def PassCode(self):
        self.Candidate['PassCode'] = 'NormalSNV'
        SubCandidate = self.Candidate.loc[self.Candidate['Genotype'].isin(['0/1', '1/1'])]
        for I in SubCandidate.index:
            currentPos = SubCandidate.loc[I, 'Pos']
            Start = currentPos - 10
            End = currentPos + 10
            SubDf = SubCandidate.loc[(SubCandidate['Pos']>=Start)&(SubCandidate['Pos']<=End)]
            altBase = SubCandidate.loc[I, 'Alt']
            if altBase not in ['A', 'T', 'C', 'G']:
                SubCandidate.loc[I, 'PassCode'] = 'Indel'
            if SubDf.shape[0] > 1:
                if self.Candidate.loc[I, 'PassCode'] == 'NormalSNV':
                    self.Candidate.loc[I, 'PassCode'] = 'ClusterSNV'
    
    def BaseCaller(self):
        '''
        genotyping for each point
        '''
        # Label SNV candidate with too many Error 
        self.FillError()
        # Cut the result dataframe to ensure all Candidate have neighborhood information
        self.Candidate = self.Candidate.loc[(self.Candidate['Pos']>=self.Start)&(self.Candidate['Pos']<=self.End)]      
        self.ABdf = self.ABdf.loc[(self.ABdf['Pos']>=self.Start)&(self.ABdf['Pos']<=self.End)]
        JudgeDf = self.Candidate[['Ref','RefCount','Alt','AltCount', 'Decision','DecisionGT','AlleleDropOut']]
        JudgeDf['Genotype'] = '0/0'
        JudgeDf['BaseCall'] = ''
        # Basecalling
        JudgeDf.loc[JudgeDf['DecisionGT']=='Heterozygous', 'Genotype'] = '0/1'
        JudgeDf.loc[JudgeDf['DecisionGT']=='Heterozygous', 'BaseCall'] = JudgeDf.loc[JudgeDf['DecisionGT']=='Heterozygous', 'Ref'] + "/" + JudgeDf.loc[JudgeDf['DecisionGT']=='Heterozygous', 'Alt']
        # Homozygous SNP 0/0 or 1/1
        JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']>JudgeDf['RefCount']), 'Genotype'] = '1/1'
        JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']>JudgeDf['RefCount']), 'BaseCall'] = JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']>JudgeDf['RefCount']), 'Alt'] + "/" + JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']>JudgeDf['RefCount']), 'Alt']
        JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']<JudgeDf['RefCount']), 'BaseCall'] = JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']<JudgeDf['RefCount']), 'Ref'] + "/" + JudgeDf.loc[(JudgeDf['DecisionGT']=='Homozygous') & (JudgeDf['AltCount']<JudgeDf['RefCount']), 'Ref']
        # Select SNV position to loose the demention
        JudgeDf_NonSNV = JudgeDf.loc[(JudgeDf['Genotype']=='0/0')]     # ref/ref base 
        JudgeDf_SNVs = JudgeDf.loc[(JudgeDf['Genotype']!='0/0')]       # alt is SNV 
        ConfuseIndex = list(JudgeDf.loc[(JudgeDf['Decision']!='Homozygous')&(JudgeDf['Decision']!='Heterozygous')&(JudgeDf['Genotype']!='0/0')].index)  # alt is SNV but high pcr risk
        self.JUDGEMENT = {"Non_SNV":JudgeDf_NonSNV,
                          "SCSNV":JudgeDf_SNVs,
                          "RawJudgeDf": JudgeDf.loc[(JudgeDf['Genotype']!="0/0") | ((JudgeDf['Decision']!='Homozygous')&(JudgeDf['Decision']!='Heterozygous'))], 
                          "ConfuseIndex": ConfuseIndex}
        self.Candidate['Genotype'] = JudgeDf['Genotype']
        self.Candidate['BaseCall'] = JudgeDf['BaseCall']
        # label Indel and 10bp continuous SNVs
        self.PassCode()
        # Cut the output data (make sure all point at the center of data)
