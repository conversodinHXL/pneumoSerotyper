#!/usr/bin/env python

__author__ = 'Chrispin Chaguza'
__copyright__ = "Copyright 2015"
__license__ = "GPL"
__version__ = "1.0.0"
__email__ = "Chrispin.Chaguza@liv.ac.uk"

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#    Copyright 2014 University of Liverpool All rights reserved.
#    The views and conclusions contained in the software and documentation are those of the
#    authors and should not be interpreted as representing official policies, either expressed
#    or implied, of University of Liverpool.

#    pneumoSerotyper.py version 1.0.0
#    LAST UPDATE: January, 2015.

#    Please report any bugs to Chrispin.Chaguza@liv.ac.uk.


import os
import sys
import glob
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML
import argparse
import math
import time
import subprocess
import timeit
import threading
import shutil
import operator


programUsage = """NAME
		pneumoSerotyper.py - Tool for identifying pneumococcal serotypes with Nucmer.
SYNOPSIS
		pneumoSerotyper.py -i seq1,seq2,...,seqN -o output -c minMatchLength
DESCRIPTION
		Run Nucmer on the provided whole genome sequences and process the matches to determine
		how they similar they are to the pneumococcal capsule polysaccharise synthesis (cps) locus
		sequence which encodes the outer cell polysaccharides that determines the serotypes.

		By default the input file format for the sequences is FASTA. The user is required to
		supply the input files and optionally specify other Options. The program can be run as
		follows;

		#Simplest way to run it is to provide the input sequences in fasta format (default).
		pneumoSerotyper.py -i *.fasta -o output

		#You can optionally specify other options minimum match length to be considered
		pneumoSerotyper.py -i *.fasta -o output -c 500

		A summary of the output file names is given at the end of program's execution.
OPTIONS
		-h
		Help. Print this help information and exit.

		-i filenames
		Specify the input fasta sequence files for the program. You can use wildcard (*)
		to speficy multiple files.

		-o filename
		Specify the prefix for the output files. By default the output files will have
		a prefix "Serotypes.Summary" if not specified by the user.

		-c value
		Specify the minimum length of the matches to be considered (default:500). Note that
		smaller values will result in more spurious matches.

		-r
		Remove Nucmer output files (default is to keep the files)
AUTHOR
		Chrispin Chaguza, Chrispin.Chaguza@liverpool.ac.uk. June 2015
FILES
		pneumoSerotyper.py
DEPENDENCIES
		MUMmer: http://mummer.sourceforge.net/
		biopython: http://biopython.org/wiki/Main_Page

		If you have homebrew installed on your computer, you can install these dependencies
		using 'brew install mummer' and 'brew install biopython'.
"""

MATCH_LENGTH=500

def readUserArguments(UserArgs):
    """
	subroutine to get input files and user options
	"""

    Options = argparse.ArgumentParser(UserArgs[0],
                        description="pneumoSerotyper.py - Tool for identifying pneumococcal serotypes with Nucmer",
                        prefix_chars="-",
                        add_help=False,
                        epilog="Chrispin Chaguza (Chrispin.Chaguza@liv.ac.uk)")

    Options.add_argument("-i",
                         action="store",
                         nargs="*",
                         required=False,
                         metavar="Input",
                         dest="Input",
                         help="input genomes (fasta format)")

    Options.add_argument("-o",
                         action="store",
                         nargs=1,
                         required=False,
                         metavar="Output",
                         dest="Output",
                         help="output prefix (default=Serotypes.Summary)",
                         default="Output")

    Options.add_argument("-l",
                         action="store",
                         nargs=1,
                         required=False,
                         metavar="Match_Length",
                         dest="Match_Length",
                         help="minimum match length (default=500)",
                         default="500")

    Options.add_argument("-r",
                         action="store_true",
                         dest="Keep",
                         help="remove Nucmer output files (keep by default)")

    Options.add_argument("-h",
                         action="store_true",
                         dest="Help",
                         help="show detailed help")

    Options = Options.parse_args()

    return Options


def checkUserArguments(UserOptions):
    Options = UserOptions
    OptionsVars = {}

    if Options.Help:
        sys.stdout.write(str(programUsage) + "\n")
        sys.exit()

    if Options.Input:
        OptionsVars["i"] = Options.Input[0:]

        for i in OptionsVars["i"]:
            if not os.path.exists(i):
                showErrorMessage("problem with input files, file "+str(i)+" not found.")
                sys.exit()
    else:
        showErrorMessage("input files (-i) are required")
        sys.exit()

    if Options.Output != "Serotypes.Summary":
        OptionsVars["o"] = Options.Output[0:][0]
    else:
        OptionsVars["o"] = Options.Output[0:]

    if Options.Match_Length != "500":
        if (int(Options.Match_Length[0:][0]) > 0 ):
            OptionsVars["l"] = Options.Match_Length[0:][0]
            global MATCH_LENGTH
            MATCH_LENGTH = OptionsVars["l"]
        else:
            showErrorMessage("sequence coverage (-l) should be >0 (default=500)")
    else:
        OptionsVars["l"] = Options.Match_Length[0:]

    OptionsVars["r"] = Options.Keep

    return OptionsVars


def showErrorMessage(ErrorMessage):
    """
    display error message on the terminal and quit
    """
    sys.stdout.write("\nerror: " + str(ErrorMessage) + "\n")
    sys.stdout.write("\nuse -h option to see more detailed help\n")
    sys.exit()


def showProgramStatus(ItemList, ItemPos):
    """
    display percentage progress for the the program
    """
    NumElements = len(ItemList)
    ProgressChars = "="

    HashChars = ProgressChars * int(round(float(ItemPos + 1) / float(NumElements) * 100))
    SpaceChars = " " * int(round(100 - len(HashChars)))
    PercentChars = float(ItemPos + 1) / NumElements * 100

    if (ItemPos + 1) >= len(ItemList) or PercentChars >= 100 or \
            (int(float(ItemPos + 1) / float(NumElements) * 100) >= 100):
        sys.stdout.write("\r|{0}| {1:.2f}% {2}".format(ProgressChars * 100, 100, ""))
        sys.stdout.flush()
        sys.stdout.write("\r|{0}| {1:.2f}% ".format(ProgressChars * 100, 100))

    else:
        sys.stdout.write("\r|{0}| {1:.2f}%  {2}".format(HashChars + SpaceChars, PercentChars, ""))
        sys.stdout.flush()
        sys.stdout.write("\r|{0}| {1:.2f}%  {2}".format(HashChars + SpaceChars, PercentChars, ""))


def RunNucmerThread(nucmerCommand, nucmerOutput):
    """
    Compare each genome against reference Cps locus sequences
    """
    try:
        FNULL = open(nucmerOutput, "wb")

        nucCommand = []

        nucCommand.append("nucmer")
        nucCommand.append("-p")
        nucCommand.append(nucmerCommand[2])
        nucCommand.append(nucmerCommand[0])
        nucCommand.append(nucmerCommand[1])

        subprocess.check_call(nucCommand,stdout=FNULL,stderr=FNULL)

        deltaFilterCommand = []

        with open(nucmerCommand[2]+".filter.delta","w") as fout:
            deltaFilterCommand.append("delta-filter")
            deltaFilterCommand.append("-i")
            deltaFilterCommand.append("90")
            deltaFilterCommand.append("-l")
            deltaFilterCommand.append(str(MATCH_LENGTH))
            deltaFilterCommand.append("-q")
            deltaFilterCommand.append("-r")
            deltaFilterCommand.append(nucmerCommand[2]+".delta")

            subprocess.check_call(deltaFilterCommand,stdout=fout,stderr=FNULL)

        coordsFilterCommand = []

        with open(nucmerCommand[2]+".coords","w") as fout:
            coordsFilterCommand.append("show-coords")
            coordsFilterCommand.append("-d")
            coordsFilterCommand.append("-T")
            coordsFilterCommand.append("-l")
            coordsFilterCommand.append("-o")
            coordsFilterCommand.append("-r")
            coordsFilterCommand.append(nucmerCommand[2]+".filter.delta")

            subprocess.check_call(coordsFilterCommand,stdout=fout,stderr=FNULL)

        FNULL.close()

    except (StandardError, KeyboardInterrupt, SystemExit):
        exceptType, exceptoObj, exceptTb = sys.exc_info()

        print "\nunknown error occurred OR the user killed the program " \
              "(line #:"+str(exceptTb.tb_lineno)+")\n\n"
        sys.exit()


def main():
    """
    runs the main program and process the user supplied genomes
    """

    Args = readUserArguments(sys.argv[:])

    inputOptions = checkUserArguments(Args)

    print "---------------------------------------------------------------------------------------------------------------"
    print "-   pneumoSerotyper v1.0.0                                                                                    -"
    print "---------------------------------------------------------------------------------------------------------------"
    #print "-                                                                                                             -"
    print "-   Program for serotyping pneumococcal isolates based on assemblies or whole genome sequences (fasta)        -"
    print "-   Institute of Infection and Global Health, University of Liverpool, UK                                     -"
    print "-   All rights reserved                                                                                       -"
    #print "-                                                                                                             -"
    print "---------------------------------------------------------------------------------------------------------------"

    print "\nchecking if mummer (http://mummer.sourceforge.net/) is installed"

    checkNucmer=0
    nucmerPath=""

    for execPath in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(execPath, "nucmer")):
            nucmerPath=os.path.join(execPath, "nucmer")
            checkNucmer+=1

    if not checkNucmer:
        print "fake sure you have installed mummer/nucmer (http://mummer.sourceforge.net/)"
        sys.exit()
    else:
        print "found nucmer executable: "+nucmerPath

    serotypesData = {
        "HE651292":"NT",
        "CR931632": "1",
        "CR931633": "2",
        "CR931634": "3",
        "CR931635": "4",
        "CR931636": "5",
        "CR931637": "5",
        "CR931638": "6A",
        "CR931639": "6B",
	"HM171374": "6D",
	"EF538714": "6C",
	"KC832411": "6F",
	"KC832410": "6G",
	"CR931640": "7A",
        "CR931641": "7B",
        "CR931642": "7C",
        "CR931643": "7F",
        "CR931644": "8",
        "CR931645": "9A",
        "CR931646": "9I",
        "CR931647": "9N",
        "CR931648": "9V",
        "CR931649": "10A",
        "CR931650": "10B",
        "CR931651": "10C",
        "CR931652": "10F",
        "CR931653": "11A",
        "CR931654": "11B",
        "CR931655": "11C",
        "CR931656": "11D",
        "CR931657": "11F",
        "CR931658": "12A",
        "CR931659": "12B",
        "CR931660": "12F",
        "CR931661": "13",
        "CR931662": "14",
        "CR931663": "15A",
        "CR931664": "15B",
        "CR931665": "15C",
        "CR931666": "15F",
        "CR931667": "16A",
        "CR931668": "16F",
        "CR931669": "17A",
        "CR931670": "17F",
        "CR931671": "18A",
        "CR931672": "18B",
        "CR931673": "18C",
        "CR931674": "18F",
        "CR931675": "19A",
        "CR931676": "19B",
        "CR931677": "19C",
        "CR931678": "19F",
        "CR931679": "20",
        "CR931680": "21",
        "CR931681": "22A",
        "CR931682": "22F",
        "CR931683": "23A",
        "CR931684": "23B",
        "CR931685": "23F",
        "CR931686": "24A",
        "CR931687": "24B",
        "CR931688": "24F",
        "CR931689": "25A",
        "CR931690": "25F",
        "CR931691": "27",
        "CR931692": "28A",
        "CR931693": "28F",
        "CR931694": "29",
        "CR931695": "31",
        "CR931696": "32A",
        "CR931697": "32F",
        "CR931698": "33A",
        "CR931699": "33B",
        "CR931700": "33C",
        "CR931701": "33D",
        "CR931702": "33F",
        "CR931703": "34",
        "CR931704": "35A",
        "CR931705": "35B",
        "CR931706": "35C",
        "CR931707": "35F",
        "CR931708": "36",
        "CR931709": "37",
        "CR931710": "38",
        "CR931711": "39",
        "CR931712": "40",
        "CR931713": "41A",
        "CR931714": "41F",
        "CR931715": "42",
        "CR931716": "43",
        "CR931717": "44",
        "CR931718": "45",
        "CR931719": "46",
        "CR931720": "47A",
        "CR931721": "47F",
        "CR931722": "48"}

    checkCount=0
    checkCPSDbase=0

    print "\nchecking the pneumococcal CPS reference sequence database"

    if os.path.exists(os.path.expanduser("~") + "/CPS_DBASE/"):
        AvailableSeqs = []

        for cpsSeq in glob.glob(os.path.expanduser("~") + "/CPS_DBASE/*.fasta"):
            AvailableSeqs.append(os.path.basename(cpsSeq).split(".")[0])

        for cpsSeqPos, cpsSeq in enumerate(serotypesData.keys()):
            if os.path.basename(cpsSeq).split(".")[0] in AvailableSeqs:
                pass
            else:

                if checkCount==0:
                    print "setting up a database for the pneumococcal CPS reference sequences\n"

                print "downloading reference CPS sequence for serotype ", serotypesData[cpsSeq.split(".")[0]]

                try:
                    tmpCpsSeq = Entrez.efetch(email="your@mail.com", db="nucleotide", id=cpsSeq.split(".")[0]
                                              , retmode="text", rettype="gb")
                    cpsSeqRecord = SeqIO.read(tmpCpsSeq, "genbank")
                    SeqIO.write(cpsSeqRecord, os.path.expanduser("~") + "/CPS_DBASE/" +
                                cpsSeq.split(".")[0] + ".fasta", "fasta")

                except (StandardError, KeyboardInterrupt, SystemExit):
                    print "\n\nfailed to download reference CPS sequences (Are you connected to internet?)."
                    print "Or perhaps you cancelled the program execution.\n\n"

                    sys.exit()

                checkCount+=1
                checkCPSDbase+=1

    else:
        os.makedirs(os.path.expanduser("~") + "/CPS_DBASE/")

        for cpsSeq in serotypesData.keys():

            if checkCount==0:
                print "setting up a database for the pneumococcal CPS reference sequences\n"

            print "downloading reference CPS sequence for serotype ", serotypesData[cpsSeq.split(".")[0]]

            try:
                tmpCpsSeq = Entrez.efetch(email="your@mail.com", db="nucleotide", id=cpsSeq.split(".")[0],
                                          retmode="text", rettype="gb")
                cpsSeqRecord = SeqIO.read(tmpCpsSeq, "genbank")
                SeqIO.write(cpsSeqRecord, os.path.expanduser("~") + "/CPS_DBASE/" + cpsSeq.split(".")[0] + ".fasta",
                            "fasta")

            except (StandardError, KeyboardInterrupt, SystemExit):
                print "\nfailed to download reference CPS sequence (Are you connected to internet?)."
                print "Or perhaps you cancelled the program execution.\n\n"
                sys.exit()

            checkCount+=1
            checkCPSDbase+=1

    if checkCPSDbase==0:
        print "located pneumococcal CPS reference sequence dbase: "+os.path.expanduser("~") + "/CPS_DBASE/ [OK]"
    else:
        print "created CPS reference sequence dbase: "+os.path.expanduser("~") + "/CPS_DBASE/ [OK]"

    if not os.path.exists(inputOptions["o"] + ".Tmp.Files/"):
        os.makedirs(inputOptions["o"] + ".Tmp.Files/")
    else:
        pass

    print "\nrunning Nucmer againt the CPS reference sequences"

    nucmerOutputFiles = []
    inputGenomeNames = []

    currentThreads = threading.activeCount()

    for cpsSeqPos, cpsSeq in enumerate(serotypesData):
        inSeqFName = os.path.basename(cpsSeq).split(".")[0]

        for isolateSeqPos, isolateSeq in enumerate(inputOptions["i"]):
            isolateFName = os.path.basename(isolateSeq).split(".")[0]

            nucmerCommand = []

            nucmerCommand.append(isolateSeq)
            nucmerCommand.append(os.path.expanduser("~") + "/CPS_DBASE/"+cpsSeq+".fasta")
            nucmerCommand.append(inputOptions["o"] + ".Tmp.Files/"+inSeqFName + "." + isolateFName)

            nucmerOutputFiles.append(inputOptions["o"] + ".Tmp.Files/"+inSeqFName +
                                     "." + isolateFName + ".coords")
            inputGenomeNames.append(isolateFName)

            while True:
                if threading.activeCount() < 15:
                    try:
                        nucmerThread = threading.Thread(name=inSeqFName+"."+isolateFName,
                                        target=RunNucmerThread,args=(nucmerCommand,os.devnull))
                        nucmerThread.setDaemon(True)
                        nucmerThread.start()

                    except (StandardError, KeyboardInterrupt, SystemExit):
                        exceptType, exceptoObj, exceptTb = sys.exc_info()
                        print "\nUnknown error occurred OR the user killed the program " \
                              "(line #:"+str(exceptTb.tb_lineno)+")\n\n"
                        sys.exit()

                    break

            showProgramStatus(serotypesData,(cpsSeqPos-1)+float(isolateSeqPos)/len(inputOptions["i"]))

        showProgramStatus(serotypesData, cpsSeqPos)

    inputGenomeNames=list(set(inputGenomeNames))

    while threading.activeCount() > currentThreads:
        time.sleep(2)

    print "\n\ncalculating length of each Cps sequence"

    cpsSeqLengths = {}

    for cpsSequencePos,cpsSequence in enumerate(glob.glob(os.path.expanduser("~") + "/CPS_DBASE/*.fasta")):
            cpsSequenceRecord=SeqIO.read(open(cpsSequence,"rU"),"fasta")
            cpsSeqLengths[os.path.basename(cpsSequence).split(".")[0]]=len(cpsSequenceRecord.seq)

    print "\ndetermining Cps sequence coverage and nucleotide identity"

    serotypesOutputCSV = open(inputOptions["o"]+".Serotypes.csv","w")
    serotypesOutputTXT = open(inputOptions["o"]+".Serotypes.txt","w")

    serotypesOutputCSV.write("isolate-id")
    serotypesOutputTXT.write("isolate-id")

    for eachCpsSeq in range(1,len(serotypesData)+1):
        serotypesOutputCSV.write(",,match-"+str(eachCpsSeq)+",coverage-"+str(eachCpsSeq)+
                                 ",identity-"+str(eachCpsSeq))
        serotypesOutputTXT.write("\tmatch-"+str(eachCpsSeq)+"\tcoverage-"+str(eachCpsSeq)+
                                 "\tidentity-"+str(eachCpsSeq))

    serotypesOutputCSV.write("\n")
    serotypesOutputTXT.write("\n")

    for genomeFilePos, genomeFile in enumerate(inputGenomeNames):

        serotypeCoverage = {}
        serotypeHspIdentity = {}

        serotypesOutputCSV.write(genomeFile)
        serotypesOutputTXT.write(genomeFile)

        for nucmerFilePos, nucmerFile in enumerate(nucmerOutputFiles):

            try:
                if nucmerFile.find("."+genomeFile+".")<0:
                    continue

                nucmerRecord=[str(i).strip().split("\t") for i in open(nucmerFile,"rU")]

                tempCoverage = []
                tempIdentity = []

                cpsSerotype = serotypesData[os.path.basename(nucmerFile).split(".")[0]]
                cpsSerotypeID = os.path.basename(nucmerFile).split(".")[0]

                for nucmerRecordAlign in nucmerRecord[4:]:
                    tempCoverage.append(float(nucmerRecordAlign[5]))

                    tempIdentity.append(float(nucmerRecordAlign[6]))

                if len(tempCoverage) == 0:
                    querySeqCoverage=0
                else:
                    querySeqCoverage = (sum(tempCoverage)/float(cpsSeqLengths[cpsSerotypeID]))*100

                if len(tempIdentity) == 0:
                    matchIdentity = 0
                else:
                    matchIdentity = (sum(tempIdentity)/float(len(tempIdentity)))

                if querySeqCoverage >= 100:
                    querySeqCoverage = 100
                else:
                    querySeqCoverage=querySeqCoverage

                serotypeCoverage[str(cpsSerotype)] = querySeqCoverage
                serotypeHspIdentity[str(cpsSerotype)] = matchIdentity

            except (StandardError, KeyboardInterrupt, SystemExit):
                exceptType, exceptoObj, exceptTb = sys.exc_info()
                print "\nsomething wrong with Nucmer files OR the user killed the program " \
                      "(line #:"+str(exceptTb.tb_lineno)+")\n\n"
                sys.exit()

            showProgramStatus(inputGenomeNames,(genomeFilePos-1)+nucmerFilePos/float(len(nucmerOutputFiles)))

        sortedSerotypeCoverage = sorted(serotypeCoverage.iteritems(),key=operator.itemgetter(1),reverse=True)

        cpsID = []
        cpsCoverage = []
        cpsIdentity = []

        for serotypeID,serotypeCpsCoverage in sortedSerotypeCoverage:
            cpsID.append(str(serotypeID))
            cpsCoverage.append(str(serotypeCpsCoverage))
            cpsIdentity.append(serotypeHspIdentity[serotypeID])

        for serotypeID,serotypeCpsCoverage,serotypeCpsIdentity in zip(cpsID,cpsCoverage,cpsIdentity):
            serotypesOutputCSV.write(",,"+str(serotypeID)+","+str(serotypeCpsCoverage)+","+str(serotypeCpsIdentity))
            serotypesOutputTXT.write("\t"+str(serotypeID)+"\t"+str(serotypeCpsCoverage)+"\t"+str(serotypeCpsIdentity))

        serotypesOutputCSV.write("\n")
        serotypesOutputTXT.write("\n")

        showProgramStatus(inputGenomeNames, genomeFilePos)

    serotypesOutputCSV.close()
    serotypesOutputTXT.close()

    if inputOptions["r"]:
        print "\n\nremoving temporary Nucmer output directory: "+inputOptions["o"] + ".Tmp.Files/"
        shutil.rmtree(inputOptions["o"] + ".Tmp.Files/")
    else:
        print "\n"

    print "CPS reference sequence database: ",os.path.expanduser("~") + "/CPS_DBASE/"
    print "inferred serotypes: ",inputOptions["o"]+".Serotypes.csv"

    if not inputOptions["r"]:
        print "output file directory: ",inputOptions["o"] + ".Tmp.Files"

    print "\n\nyou can check similarities between various pneumococcal CPS loci using link below: "
    print "http://spneumoniae.mlst.net/CPS/database/v3/ACT"

    print "\n--------------------------------------------------------------------------------------------------------------"
    print "-        Finished                                                                                            -"
    print "--------------------------------------------------------------------------------------------------------------"


if __name__ == "__main__":
    """
    start the main program
    """
    main()
