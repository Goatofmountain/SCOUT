# NBJUDGE_wholeGenome main
import sys, getopt
import time
import os
import multiprocessing
from multiprocessing import Manager,Pool

# fix source path
Bin_dir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
source_calculate = os.path.join(os.path.abspath(os.path.dirname(Bin_dir)), 'source/Calculate')
source_plot = os.path.join(os.path.abspath(os.path.dirname(Bin_dir)), 'source/Plot')
stat_dir = os.path.join(os.path.abspath(os.path.dirname(Bin_dir)), 'stat')
sys.path.append(source_calculate)
sys.path.append(source_plot)
sys.path.append(stat_dir)
from Candidate import *
HelpInfo = 'python SCOUT_WholeGenome.py -N <ProjectName> -r <reference> -i <inputBam> -o <outputdir> -W <WindowSize, default=30000> -M <decay mu value, default=ln(3)/10> -P <Number of CPU> -c <AimChrom> -m <runningModel>'

def WorkPip(refFile, samfile,CellName, WindowSize,mu, chrom, Start, End):
    '''
    Run Estimate pipline on a specific region
    Return the Whole object PipEST
    '''
    region ='%s:%s-%s' % (chrom, Start, End)
    print('%s: Region %s Start to Estimate' % (time.ctime(), region))
    PipEST = EstModel(bamFile=samfile,
                      ref=refFile,
                      CellName=CellName,
                      WindowSize=WindowSize,
                      mu=mu)
    PipEST.MakeCandidateDf(chrom, Start, End)
    if PipEST.Candidate.shape[0] > 0:
        PipEST.EstimateError()
        PipEST.BaseCaller()
    Summary = {"Candidate":PipEST.Candidate, 'ABdf':PipEST.ABdf}
    return(Summary)

def RunPip(Parameter):
    # Run pipeline
    ResultPool = {}
    WorkPool = Pool(processes=Parameter['Process'])
    RegionLen = Parameter['End'] - Parameter['Start']
    tmpStart = Parameter['Start']
    tmpEnd = Parameter['Start'] + 2000000
    while Parameter['End'] - tmpStart >= 2000000:
        regionName = '%s,%s,%s' % (Parameter['chrom'], tmpStart, tmpEnd)
        ResultPool[regionName] = WorkPool.apply_async(WorkPip,
                                                      args=(Parameter['refFile'],
                                                            Parameter['samfile'],
                                                            Parameter["PName"],
                                                            Parameter['WindowSize'],
                                                            Parameter['mu'],
                                                            Parameter['chrom'],
                                                            tmpStart,
                                                            tmpEnd,)
                                                     )
        tmpStart = tmpEnd + 1
        tmpEnd = tmpStart + 2000000
    if tmpStart < Parameter['End']:
        regionName = '%s,%s,%s' % (Parameter['chrom'], tmpStart, Parameter['End'])
        ResultPool[regionName] = WorkPool.apply_async(WorkPip,
                                                      args=(Parameter['refFile'],
                                                            Parameter['samfile'],
                                                            Parameter["PName"],
                                                            Parameter['WindowSize'],
                                                            Parameter['mu'],
                                                            Parameter['chrom'],
                                                            tmpStart,
                                                            Parameter['End'])
                                                     )
    WorkPool.close()
    WorkPool.join()
    CandidateList = []
    ABdfList = []
    for k in sorted(list(ResultPool.keys()), key=lambda s:int(s.split(",")[1])):
        print('%s: Start to fetch region : %s' % (time.ctime(), k))
        res = ResultPool[k].get()
        if not res['Candidate'].empty:
            CandidateList.append(res['Candidate'])
            ABdfList.append(res['ABdf'])
    Candidate = pd.concat(CandidateList, axis=0)
    ABdf = pd.concat(ABdfList, axis=0)
    os.chdir(Parameter['outputdir'])
    print('%s: Now Write out the Result' % time.ctime())
    Candidate.to_csv("%s.%s.CandidateDf.csv" % (Parameter['PName'], Parameter['chrom']))
    ABdf.to_csv('%s.%s.ABdf.csv' % (Parameter['PName'], Parameter['chrom']))
    print('%s: %s.%s Finished' % (time.ctime(), Parameter['PName'], Parameter['chrom']))
    # Write out the vcf File 
    CandidateDf = Candidate
    VCFFile = "%s.%s.vcf" % (Parameter['PName'], Parameter['chrom'])
    CellName = Parameter['PName']
    bamFile = Parameter['samfile']
    refFile = Parameter['refFile']
    CallVCF(CandidateDf, VCFFile, refFile, bamFile, CellName)
    print("%s: VCF caller for %s.%s Finished!" % (time.ctime(), Parameter['PName'], Parameter['chrom']))
    del WorkPool
    
def main(argv):
    Parameter = {"samfile":"",
                 "outputdir":os.getcwd(),
                 "chrom":'chr1',
                 "Start":100000000,
                 "End":103000000,
                 "WindowSize":30000,
                 "mu":np.log(3)/10,
                 "PName":None,
                 "refFile":None,
                 "Process":1, 
                 "model":"SpecificRegion"
                }
    TableFile = None
    RequireParm = ['samfile', 'chrom', 'Start', 'End', 'PName', 'refFile', 'Process']
    try:
        opts,args = getopt.getopt(argv, "hN:r:i:o:c:S:E:W:M:P:m:", ['Name=', 'ref=', 'ifile=', 'ofile=', 'chr=', 'Start=', 'End=', 'WindowSize=', 'Mu=', "Process=", "model="])
    except getopt.GetoptError:
        print(HelpInfo)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(HelpInfo)
            sys.exit(2)
        elif opt in ('-N', '--Name'):
            Parameter['PName'] = arg
        elif opt in ('-r', '--ref'):
            Parameter['refFile'] = arg
        elif opt in ("-i", "--ifile"):
            Parameter['samfile'] = arg
        elif opt in ("-o", "--ofile"):
            Parameter['outputdir'] = arg
        elif opt in ("-c", "--chr"):
            Parameter['chrom'] = arg
        elif opt in ("-S", '--Start'):
            Parameter['Start'] = int(arg)
        elif opt in ("-E", '--End'):
            Parameter['End'] = int(arg)
        elif opt in ("-W", '--WindowSize'):
            Parameter['WindowSize'] = int(arg)
        elif opt in ("-M", '--Mu'):
            Parameter['mu'] = float(arg)
        elif opt in ("-P", "--Process"):
            Parameter['Process'] = int(arg)
        elif opt in ('-m', '--model'):
            Parameter['model'] = arg
    # Check required Parameters
    for P in RequireParm:
        if not Parameter[P]:
            print("The parametar %s is not defined please set it !" % P)
            sys.exit(2)
    Parameter['Start'] = int(Parameter['Start'])
    Parameter['End'] = int(Parameter['End'])
    Parameter['WindowSize'] = int(Parameter['WindowSize'])
    Parameter['Process'] = int(Parameter['Process'])
    PName = Parameter['PName']
    print('%s: %s Start' % (time.ctime(), PName))
    resultDir = os.path.join(Parameter['outputdir'], Parameter['PName'])
    if not os.path.exists(resultDir):
        os.mkdir(resultDir)
    Parameter['outputdir'] = resultDir
    if Parameter['model'] == 'WholeGenome':
        chromList = Parameter['chrom'].split(",")
        try:
            f = open(Parameter['refFile'] + ".fai", 'r')
        except:
            print("reference file do not have fai index, please build an index using samtools faidx")
            sys.exit(2)
        with open(Parameter['refFile'] + ".fai") as f:
            chrLen = {}
            for records in f.readlines():
                chrLen[records.strip().split("\t")[0]] = int(records.strip().split("\t")[1])
        for AimChrom in chromList:
            if AimChrom in chrLen.keys():
                End = chrLen[AimChrom]
                Start = 1
                tmpParameter = Parameter
                tmpParameter['Start'] = Start
                tmpParameter['End'] = End
                tmpParameter['chrom'] = AimChrom
                RunPip(tmpParameter)
            elif AimChrom not in chrLen.keys():
                print('%s not in the chromosome list. Skipping ...' % AimChrom)
    else:
        RunPip(Parameter)

if __name__ == '__main__':
    main(sys.argv[1:])
