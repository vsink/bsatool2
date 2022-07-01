package main

import (
	"bsatoolMod/src/bsatcalc"
	"bsatoolMod/src/bsatdb"
	"bsatoolMod/src/bsatf"
	"bsatoolMod/src/bsatseq"
	"github.com/pterm/pterm"
	"gopkg.in/alecthomas/kingpin.v2"

	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bsatoolMod/src/bsatvcf"
	// "bufio"

	"runtime/pprof"
	// "time"

	"fmt"

	"log"

	// "encoding/gob"

	// "net/http"
	"os"

	// "reflect"
	// "regexp"

	"strings"
	// "./amino"
	// "./amino"
	// "html/template"
	// "github.com/pkg/browser"
)

const (
// 	logo = `
//  _______    ________     __  ___________  ______      ______    ___
// |   _  "\  /"       )   /""\("     _   ")/    " \    /    " \  |"  |
// (. |_)  :)(:   \___/   /    \)__/  \\__/// ____  \  // ____  \ ||  |
// |:     \/  \___  \    /' /\  \  \\_ /  /  /    ) :)/  /    ) :)|:  |
// (|  _  \\   __/  \\  //  __'  \ |.  | (: (____/ //(: (____/ //  \  |___
// |: |_)  :) /" \   :)/   /  \\  \\:  |  \        /  \        /  ( \_|:  \
// (_______/ (_______/(___/    \___)\__|   \"_____/    \"_____/    \_______)
// ` +
// 		`BSATool - Bacterial Snp Annotation Tool ` + "\n" +
// 		`      Laboratory of Social and Epidemic Infections
//  Scientific Centre for Family Health and Human Reproduction Problems
//      	(c) V.Sinkov, Irkutsk, Russia, 2017-2022
//
// 	`
// list   = "list"
// ncFlag = "NC"
// aaFlag = "AA"
)

var (

	//
	// Database flags
	version string
	build   string

	// app       = kingpin.New(logo, "BSATool - Bacterial Snp Annotation Tool")
	// appAuthor = app.Author("V.Sinkov")
	// appVer    = app.Version(fmt.Sprintf("%v %v", version, build))

	gbWeb     = kingpin.Flag("web", " Open results in web browser").Short('w').Bool()
	gbXLSX    = kingpin.Flag("xlsx", " Export to XLSX format").Short('x').String()
	gbVerbose = kingpin.Flag("verbose", "Show additional information ").Short('v').Default("false").Bool()
	gbIndex   = kingpin.Flag("index", "Calculate Complex Index for amino acid changes").Default("false").Bool()
	gbPort    = kingpin.Flag("port", "Use your own localhost:port (default:8080)").Default("8080").String()
	gbNoSeq   = kingpin.Flag("noseq", "Don't show nucleotides").Default("false").Bool()
	gbDebug   = kingpin.Flag("debug", "Debug mode").Default("false").Bool()
	// gbLog            = kingpin.Flag("log", "write log file").Default("false").Bool()
	gbExcludeGenes   = kingpin.Flag("exclude-genes", "file with genes which should be excluded from mkseq").String()
	gbExcludeRegions = kingpin.Flag("exclude-regions", "file with coordinates(start end)  which should be excluded from mkseq").String()
	gbExcludeSnp     = kingpin.Flag("exclude-snp", "file with genes which should be excluded from mkseq").String()
	gbRandomize      = kingpin.Flag("randomize", "Set on randomizer").Default("false").Bool()
	mkdb             = kingpin.Command("mkdb", "Create database")
	dbName           = mkdb.Flag("out", "Name of database").Short('o').Required().String()
	dbGenbank        = mkdb.Flag("gb", "Name of genbank file").Short('i').Required().String()
	gbDP             = kingpin.Flag("dp", "dp value for filtering vcfs").Default("1").Int()
	gbIGR            = kingpin.Flag("igr", "Show igr regions in results").Default("false").Bool()
	gbInDel          = kingpin.Flag("indel", "indel detection").Bool()
	// gbAbout        = kingpin.Flag("about", "About programm").Bool()

	// Annotation flags

	annAction = kingpin.Command("annotate", "This is command to allow annotate VCF file").Alias("annotation").Alias("ann").Alias("a")
	annDB     = annAction.Flag("db", "Name of database").Short('b').Required().String()
	annVCF    = annAction.Flag("vcf", "Input VCF file").Short('i').Required().String()
	// annWeb           = annAction.Flag("web", "").Bool()
	annMakeSeq       = annAction.Flag("mkseq", "NC, AA, Binary").Short('m').String()
	annOutFormat     = annAction.Flag("output-format", "Output format for binary data (phylip, nexus)").String()
	annMakeSeqRef    = annAction.Flag("ref", "Generate reference sequence").Short('r').Default("false").Bool()
	annWithFilenames = annAction.Flag("wfn", "Show filenames in list annotated VCF's").Short('n').Bool()
	annBench         = annAction.Flag("annprof", "cpuprofile").String()
	annSeqLen        = annAction.Flag("len", "length of sequence").Int()
	annShowFileName  = annAction.Flag("filename", "Show filename in first column for next analysis").Default("false").Bool()
	annPosFile       = annAction.Flag("pos-file", "Create positions.txt file with information of each pos in sequence").String()
	annMinPos        = annAction.Flag("min-pos", "Set minimal threshhold of presence of position in sequence").Default("1").Int()
	annMaxPos        = annAction.Flag("max-pos", "Set maximal threshhold of presence of position in sequence").Int()
	// annAltPosPercent     =annShowFileName annAction.Flag("alt-percent", "Maximal percent of presence of reference allele in sequence (99 default)").Default("1").Int()
	annAltRange = annAction.Flag("range", "Min-Max range values (percent) of presence alternative nucleotide in sequence").Default("10:90").String()
	// annBench         = annAction.Flag("cpuprofile", "cpuprofile").String()

	// compute statistic options

	statAction = kingpin.Command("stat", "Calculates statistic tests").Alias("s")
	statDB     = statAction.Flag("db", "Database file").Short('b').Required().String()
	statTask   = statAction.Flag("action", "Type of action:share, snp,dnds, ginfo,bed, matrix,range, circos").Short('a').Required().String()
	// statWeb     = annAction.Flag("web", "").Short('w').Bool()
	statInFile     = statAction.Flag("in", "Input file").Short('i').String()
	statInFileList = statAction.Flag("flist", "Input file").Short('f').String()
	statOutFile    = statAction.Flag("out", "Output file").Short('o').String()
	statTypeOf     = statAction.Flag("type", "Type of matrix (binary, gc3, dnds, nc, locus. freq, jw, summary").Short('t').String()
	statInRule     = statAction.Flag("rule", "Input rule file").Short('r').String()
	statAll        = statAction.Flag("all", "show all dNdS results").Default("false").Bool()
	statBench      = statAction.Flag("statprof", "cpuprofile").String()
	statMakeSeq    = statAction.Flag("mkseq", "make seq from snp filelist").Bool()
	// statRefGenomeName   = statAction.Flag("ref", "set custom genome name").String()
	statCircosTypeOf    = statAction.Flag("typeof", "make circos file without IGR regions").String()
	statBedTypeOf       = statAction.Flag("bed-typeof", "cds, igr").String()
	statCircosBandColor = statAction.Flag("color", "set color to band in circos").String()
	statGenomeName      = statAction.Flag("genome-name", "genome name").String()
	statShowAnnotation  = statAction.Flag("annotation", "show annotations to genes").Bool()
	statNbrOfSNP        = statAction.Flag("snp-number", "number of snp for MST matrix").Int()
	statMinLocusCount   = statAction.Flag("min-locus-count", "the minimal number of snp per locus ").Int()
	statGroupFromFile   = statAction.Flag("group", "File with filenames and their groups").String()
	statVCF             = statAction.Flag("vcf", "Input VCF file").String()
	statLocus           = statAction.Flag("locus", "locus name").String()
	statDnDsCountTh     = statAction.Flag("th", "DnDs threshold (How many files consists nonzero values (Default > 1))").Int()
	// statGStrain         = statAction.Flag("strain", "set strain name").String()
	statMkSeq    = statAction.Flag("with_seq", "make sequence").Bool()
	statCircos   = statAction.Flag("circos", "export to circos file").Bool()
	statCircosGC = statAction.Flag("gc", "export to circos file").Bool()
	// statGName           = statAction.Flag("gname", "set strain name").String()
	statFlankLeft  = statAction.Flag("flank_left", "number of left nucleotides").Int()
	statFlankRight = statAction.Flag("flank_right", "number of right nucleotides").Int()
	statAsTable    = statAction.Flag("table", "show results of check info as table").Bool()
	statAsGene     = statAction.Flag("gene", "show mutation as gene coordinates").Bool()
	statReverse    = statAction.Flag("reverse", "make seq for complementnary genes in gene direction").Bool()
	// statAsBinary        = statAction.Flag("binary", "show results of check info as binary sequence").Bool()
	statTH       = statAction.Flag("nbr_pos", "nbr found SNP (default = nbr files-1)").Int()
	statDnDsStat = statAction.Flag("show-info", "show total statistics per genome").Bool()
	// statAsNuc = statAction.Flag("nuc", "make table with check range command as nuc").Bool()

	// statMkSeq   = statAction.Flag("mkseq", "").Bool()
	// statTH      = statAction.Flag("th", "").Int()

	filterAction = kingpin.Command("filter", "FilterVcfFILEs")
	dpMapAction  = kingpin.Command("dpmap", "DP mapping")
	dpMapDB      = dpMapAction.Flag("db", "Name of database").Short('b').Required().String()
	// filterDB     = filterAction.Flag("db", "Database file [Created by mkdb command]").Short('b').Required().String()

	infoAction     = kingpin.Command("info", "Get information")
	infoLocus      = infoAction.Flag("locus", "locus name").String()
	infoDB         = infoAction.Flag("db", "Database file [Created by mkdb command]").Short('b').Required().String()
	infoRanges     = infoAction.Flag("range", "Show sequence of nucleotides [start:end]").String()
	infoCodons     = infoAction.Flag("codons", "Show sequence of codons [start:end]").String()
	infoShowAs     = infoAction.Flag("showas", "Show as:direct (from left to right),gene(direction as in gene) ").String()
	infoMask       = infoAction.Flag("mask", "pos:ref:alt").String()
	infoReplace    = infoAction.Flag("replace", "pos:ref:alt").String()
	infoFlankLeft  = infoAction.Flag("flank_left", "number of left nuclotides").Int()
	infoFlankRight = infoAction.Flag("flank_right", "number of right nuclotides").Int()
	infoIUPAc      = infoAction.Flag("iupac", "pos:ref:alt").Bool()

	devAction = kingpin.Command("dev", "Developer mode.").Alias("debug")
	devDB     = devAction.Flag("db", "Database file").Short('b').Required().String()
	devTask   = devAction.Flag("action", "Action...").Short('a').Required().String()
	devPwd    = devAction.Flag("pwd", "Password").Short('p').Required().String()

	betaAction = kingpin.Command("beta", "Beta test mode for testing new functions")
	betaTask   = betaAction.Flag("action", "Action...").Short('a').Required().String()
	betaDB     = betaAction.Flag("db", "Database").Short('b').Required().String()
	betaInFile = betaAction.Flag("in", "Input file").Short('i').String()
	// betaOutFile = betaAction.Flag("out", "Output file").Short('o').String()

	pileup2fasta  = kingpin.Command("parse_pileup", "Pileup parser")
	pileupInFile  = pileup2fasta.Flag("in", "Input file").Short('i').String()
	pileupTH      = pileup2fasta.Flag("th", "DnDs threshold (How many files consists nonzero values (Default > 1))").Int()
	pileupGStrain = pileup2fasta.Flag("strain", "set strain name").String()
	pileupMkSeq   = pileup2fasta.Flag("mkseq", "make sequence").Bool()
	pileupShowSeq = pileup2fasta.Flag("showseq", "make sequence").Bool()
	pileupGName   = pileup2fasta.Flag("genome", "set genome name").String()
	pileupCircos  = pileup2fasta.Flag("circos", "show as circos file format").Bool()
	pileupDB      = pileup2fasta.Flag("db", "Database file").Short('b').Required().String()
	pileupMinLen  = pileup2fasta.Flag("min-len", "Минимальная длина фрагмента").Int()
)

func main() {
	// загружаем список vcf файлов
	bsatstruct.ListOfFiles = bsatf.GetListVcf()
	pterm.DefaultCenter.Println("\n")
	s, _ := pterm.DefaultBigText.WithLetters(pterm.NewLettersFromString("BSATool")).Srender()
	pterm.DefaultCenter.Println(s) // Print BigLetters with the default CenterPrinter
	pterm.DefaultCenter.WithCenterEachLineSeparately().Println(
		"-= BSATool - Bacterial Snp Annotation Tool =-\nThe Laboratory of Social and Epidemic Infections\nThe Group of Genomic Research and"+
			" Bioinformatics\nScientific Centre for Family"+
			" Health"+
			" and Human"+
			" Reproduction"+
			" Problems\n(c) V.Sinkov, Irkutsk, Russia, 2017-2022\nVersion:", version)

	// парсинг флагов
	defer os.Exit(0)

	// flag.Parse()
	// kingpin.New(logo, "BSATool - Bacterial Snp Annotation Tool")
	kingpin.Version(fmt.Sprintf("%v %v", version, build)).Author("V.Sinkov")
	kingpin.UsageTemplate(kingpin.LongHelpTemplate)

	kingpin.Parse()
	bsatstruct.Flag = bsatstruct.Flags{AnnSeqLen: *annSeqLen, AnnMakeSeqRef: *annMakeSeqRef, AnnBench: *annBench, AnnMakeSeq: *annMakeSeq,
		AnnShowFileName: *annShowFileName, AnnDB: *annDB, AnnPosFile: *annPosFile, AnnAltRange: *annAltRange, AnnMinPos: *annMinPos,
		AnnOutFormat: *annOutFormat, AnnVCF: *annVCF, AnnWithFilenames: *annWithFilenames, BetaDB: *betaDB, BetaInFile: *betaInFile, BetaTask: *betaTask,
		GbIGR: *gbIGR, GbRandomize: *gbRandomize, GbIndex: *gbIndex, GbInDel: *gbInDel, GbDP: *gbDP, GbPort: *gbPort, GbDebug: *gbDebug, GbVerbose: *gbVerbose,
		GbExcludeGenes: *gbExcludeGenes, GbExcludeRegions: *gbExcludeRegions, GbExcludeSnp: *gbExcludeRegions, GbNoSeq: *gbNoSeq, GbWeb: *gbWeb, GbXLSX: *gbXLSX,
		DbGenbank: *dbGenbank, DbName: *dbName, DevDB: *devDB, InfoCodons: *infoCodons, InfoDB: *infoDB, InfoFlankLeft: *infoFlankLeft,
		InfoFlankRight: *infoFlankRight, InfoIUPAc: *infoIUPAc, InfoLocus: *infoLocus, InfoMask: *infoMask, InfoRanges: *infoRanges, InfoShowAs: *infoShowAs,
		PileupCircos: *pileupCircos, PileupDB: *pileupDB, PileupGName: *pileupGName, PileupGStrain: *pileupGStrain, PileupInFile: *pileupInFile,
		PileupTH: *pileupTH, PileupMkSeq: *pileupMkSeq, PileupShowSeq: *pileupShowSeq, PileupMinLen: *pileupMinLen, StatFlankRight: *statFlankRight,
		StatFlankLeft: *statFlankLeft, StatTH: *statTH, StatReverse: *statReverse, StatTask: *statTask, StatTypeOf: *statBedTypeOf, StatAsTable: *statAsTable,
		StatCircosTypeOf: *statCircosTypeOf, StatBedTypeOf: *statBedTypeOf, StatCircos: *statCircos, StatDnDsCountTh: *statDnDsCountTh, StatInRule: *statInRule,
		StatAsGene: *statAsGene, StatDB: *statDB, StatDnDsStat: *statDnDsStat, StatAll: *statAll, StatBench: *statBench, StatCircosGC: *statCircosGC,
		StatCircosBandColor: *statCircosBandColor, StatGenomeName: *statGenomeName, StatInFile: *statInFile, StatInFileList: *statInFileList,
		StatGroupFromFile: *statGroupFromFile, StatLocus: *statLocus, StatOutFile: *statOutFile, StatMakeSeq: *statMakeSeq,
		StatShowAnnotation: *statShowAnnotation, StatMkSeq: *statMkSeq, StatVCF: *statVCF, StatNbrOfSNP: *statNbrOfSNP, StatMinLocusCount: *statMinLocusCount,
		DevPwd: *devPwd, DevTask: *devTask, InfoReplace: *infoReplace, AnnMaxPos: *annMaxPos, Version: version}

	switch kingpin.Parse() {

	// bsatool mkdb/create/makedb -i /--genbank FILE --out /-o FILE
	case "mkdb":

		if _, err := os.Stat(*dbGenbank); os.IsNotExist(err) {
			pterm.Error.Printf("The %v file is not exist!\n", *dbGenbank)
			os.Exit(3)
		}

		bsatdb.CreateDB(*dbGenbank, *dbName)

	case "annotate":

		// Загружаем базу данных
		bsatdb.LoadDB(*annDB)

		var (
			exGenes      = make(map[int]int)
			exGenesRaw   = make(map[int]int)
			exRegions    = make(map[int]int)
			exGenesBySNP = make(map[int]int)
			exSNPs       = make(map[int]int)
		)

		if *gbExcludeGenes != "" {
			if _, err := os.Stat(*gbExcludeGenes); os.IsNotExist(err) {
				fmt.Println("No file with excluded genes was found")
				os.Exit(3)

			} else {

				exGenesRaw = bsatf.LoadExclGenes(*gbExcludeGenes)

			}
		}

		if *gbExcludeRegions != "" {
			if _, err := os.Stat(*gbExcludeRegions); os.IsNotExist(err) {
				fmt.Println("No file with excluded regions was found")
				os.Exit(3)

			} else {

				exRegions = bsatf.LoadExclRegion(*gbExcludeRegions)
			}
		}

		if *gbExcludeSnp != "" {
			if _, err := os.Stat(*gbExcludeSnp); os.IsNotExist(err) {
				fmt.Println("No file with excluded SNPs was found")
				os.Exit(3)

			} else {

				exSNPs = bsatf.LoadExclSNP(*gbExcludeSnp)
				exGenesBySNP = bsatf.ExclGeneByPos(*gbExcludeSnp)

			}

		}
		if len(exGenesRaw) != 0 {
			for k, v := range exGenesRaw {
				exGenes[k] = v
			}
		}
		if len(exRegions) != 0 {
			for k, v := range exRegions {
				exGenes[k] = v
			}
		}
		if len(exGenesBySNP) != 0 {
			for k, v := range exGenesBySNP {
				exGenes[k] = v
			}
		}

		// fmt.Println(exGenes)

		if *annVCF == "list" || *annVCF == "*" || *annVCF == "all" {
			if *gbWeb == false && *annMakeSeq == "" {

				bsatseq.ParserBulkVCF(*annWithFilenames)
				// go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core

			} else if *gbWeb == true && strings.ToUpper(*annMakeSeq) == "NC" {
				bsatseq.NucWebServer(*gbPort, exGenes, exSNPs)
			} else if *gbWeb == false && *annMakeSeq == "NC" && strings.ToUpper(*annOutFormat) != "NEXUS" {
				// seq := MakeSeqFasta(*annMakeSeq, *gbVerbose, *annMakeSeqRef)

				/* alt-range значит, что процент встречаемости альтернативного аллеля в последовательности будет в диапазоне --alt-range=10:80
								go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core --debug --alt-range=10:80
				makeSeq
				*/

				seq := bsatseq.MakeSeqFasta(*annMakeSeq, *gbVerbose, *annMakeSeqRef, exGenes, exSNPs, *gbRandomize)
				// fmt.Println(len(seq))
				for i := 0; i < len(seq); i++ {
					fmt.Println(seq[i].Seq)
				}
				if *gbDebug == true {
					fmt.Println(seq[0].UsedPositions)
				}
			} else if *gbWeb == false && strings.ToUpper(*annMakeSeq) != "BINARY" && strings.ToUpper(*annOutFormat) == "NEXUS" {
				// fmt.Println(exGenes,exSNPs)
				bsatseq.MakeSeqNex(*annMakeSeq, *gbVerbose, exGenes, exSNPs, *gbRandomize)

			} else if *gbWeb == false && strings.ToUpper(*annMakeSeq) == "BINARY" {
				// go run bsatool.go annotate --vcf list --mkseq=NC  --db test_core
				// go run bsatool.go annotate --vcf list --mkseq=Binary  --db test_core
				// go run  bsatool.go annotate --db test_core --vcf list --exclude-genes=ppe.genes --exclude-snp=drugs3.
				// txt --mkseq=Binary --dp=30 --exclude-regions=exgenes.txt --min-pos-nbr=2 --output-format="nexus"

				bsatseq.MakeSeqBin(*annMakeSeq, *gbVerbose, exGenes, exSNPs, *gbRandomize, *gbDP)
				// bsatool annotate --db sars --vcf list  --mkseq=Alignment --indel

			} else if *gbWeb == false && strings.ToUpper(*annMakeSeq) == "ALIGNMENT" {
				bsatseq.MakeAlign(false)
			}
		} else {
			// go run bsatool.go annotate --vcf 161_RuU_m.vcf  --db test_core -w

			// qSNP := &bsatvcf.VcfQuery{File: *annVCF, OutChan: make(chan bsatvcf.VcfInfoQuery)}
			// go qSNP.Request()
			// // snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
			// snpRes := <-qSNP.OutChan
			snpsChan := make(chan []bsatstruct.SnpInfo)

			go func() {
				snpsChan <- bsatvcf.GetSNPList(*annVCF)
			}()
			snpRes := <-snpsChan

			if _, err := os.Stat(*annVCF); os.IsNotExist(err) {
				fmt.Printf("The %v file is not exist!\n", *annVCF)
				os.Exit(3)
			}
			if *gbWeb == false && *annMakeSeq == "" && *gbXLSX == "" {

				for i := 0; i < len(snpRes); i++ {
					bsatservice.PrintInConsole(snpRes[i], *gbIGR)
				}

				// parserVCF(*annVCF, true, bsatstruct.FullGenesInfo)

			} else if *gbWeb == true && *annMakeSeq == "" {
				// snps := parserVCF(*annVCF, false, bsatstruct.FullGenesInfo)
				bsatservice.PrintInWEB(snpRes, *gbPort)
			}
			if *gbXLSX != "" {
				// snps := parserVCF(*annVCF, false, bsatstruct.FullGenesInfo)
				bsatservice.ExportToExcel(snpRes, *gbXLSX)

			}
		}
	case "stat":
		cpuprofile := *statBench
		if cpuprofile != "" {
			f, err := os.Create(cpuprofile)
			if err != nil {
				log.Fatal(err)

			}
			err = pprof.StartCPUProfile(f)
			if err != nil {
				log.Fatal(err)
			}
			defer pprof.StopCPUProfile()
		}

		bsatdb.LoadDB(*statDB)

		switch *statTask {
		case "share":
			/*
			 go run bsatool.go   stat -b test_core  -a share  -v
			*/

			bsatvcf.GetShareSNP(*gbVerbose, *gbIGR, *gbWeb, bsatstruct.ListOfFiles)

		case "snp":
			/*
				go run bsatool.go   stat -b test_core  -a snp  -v
			*/
			bsatcalc.MakeSNPStatistic()
		case "dnds":

			/*
				go run bsatool.go   stat -b test_core  -a dnds  -i 161_RuU_m.vcf -v
			*/

			if _, err := os.Stat(*statInFile); os.IsNotExist(err) {
				fmt.Println("No input file found")
				os.Exit(3)
			}

			bsatcalc.CalcDnDs(*statInFile)
		case "ginfo":
			/*
				go run bsatool.go   stat -b test_core  -a ginfo
			*/
			bsatseq.TestGeneInfo(bsatstruct.FullGenesInfo)
		case "circos":
			/*
				go run bsatool.go   stat -b test_core  -a circos --typeof=[cds,igr,gene] --color=lblue --genome NC_000962
			*/
			bsatseq.ToCircos(bsatstruct.FullGenesInfo)
		case "bed":
			/*
				go run bsatool.go stat --db test_core -a bed --genome-name="NC_000962.3" --bed-typeof="CDS"
			*/
			// fmt.Println(*statGenomeName, "!!!!")
			bsatseq.ToBED(bsatstruct.FullGenesInfo, *statGenomeName)
		case "matrix":

			/*
				 go run bsatool.go stat -a matrix --db test_core -t binary (binary, gc3, dnds, table, nc, locus. freq, jw, summary) -o test.csv
				 go run bsatool.go stat   -b test_core  -a matrix -t binary  --exclude-genes=ppe.genes --exclude-snp=drugs3.txt -o test.csv --nbr_pos 2 --debug|grep "passed"
				go run bsatool.go stat  -b test_core  -a matrix -t locus -o test.csv --min-locus-count 3
			*/

			if *statTypeOf != "" && *statOutFile != "" {
				bsatcalc.MakeMatrix(*statTypeOf, *statOutFile, *gbVerbose)
			} else if *statTypeOf != "" && *statOutFile == "" {
				fmt.Println("Please, use key -o (--out) to save data.")
			}
		case "range":
			/*
				go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt
				go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt
			*/

			if *statInFile != "" {
				file := *statInFile
				res := bsatservice.GetRangeFromFile(file, *gbVerbose, *gbNoSeq, *statGenomeName)
				if *statReverse == true {
					for _, val := range res {
						fmt.Println(bsatservice.GetRevComplement(val.Seq))
					}

					bsatservice.PrintSeqRange(res, *gbWeb, *gbPort, *statCircos)
				} else {

					bsatservice.PrintSeqRange(res, *gbWeb, *gbPort, *statCircos)
				}

			}

			/*
				go run bsatool.go stat -a range --db test_core -i BWAIntersections.txt
			*/

		case "check_range":

		case "check":

			/*
				go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w
				 go run bsatool.go stat --db test_core -a check -i drugs4.txt --table
				 go run bsatool.go stat --db test_core -a check -i drugs4.txt --mkseq
			*/
			// bsatstruct.StatInRule = *statInRule

			if bsatstruct.Flag.StatMakeSeq == true {
				fmt.Println(bsatstruct.Flag.StatInFile)
				bsatf.MakeSeqFromSNPListFile(bsatstruct.Flag.StatInFile)

			} else {
				var useRule bool
				if bsatstruct.Flag.StatInRule != "" {
					useRule = true
					bsatstruct.RulesArr = bsatf.CheckRuleFromFile(bsatstruct.Flag.StatInRule)
				} else {
					useRule = false
				}
				if bsatstruct.Flag.StatInFile != "" {

					bsatf.CheckSNPfromFile(bsatstruct.Flag.StatInFile, bsatstruct.Flag.GbVerbose, bsatstruct.Flag.GbWeb, useRule)
				}
			}

		case "rule":
			/*
				go run bsatool.go stat -a check --db test_core -i drugs2.txt -r rule.txt -w
			*/

			if *statInRule != "" {
				bsatstruct.RulesArr = bsatf.CheckRuleFromFile(bsatstruct.Flag.StatInRule)

			}

		// создает fasta файл основываясь на снипах в границах указзаного диапазона
		case "Pos2Seq":
			/*  go run bsatool.go stat  -b test_core -a Pos2Seq --vcf=list -i coordsToMakeSeq.txt
			/*  go run bsatool.go stat  -b test_core -a Pos2Seq --vcf=list -i coordsToMakeSeq.txt -f CladaA

			Rv0018c	21637	23181 <-structure file
			*/

			if *statInFile != "" && *statVCF == "list" {

				if len(*statInFileList) != 0 {
					flist := bsatf.LoadFileNames(*statInFileList)
					// fmt.Println(flist)
					for _, file := range flist {
						if bsatseq.CheckVCFFormat(file) == true {
							result := bsatseq.Pos2Seq(*statInFile, file)
							locFile := strings.Replace(result.VcfFile, ".vcf", "", -1)
							fmt.Printf(">%v_%v(%v_%v)\n%v\n", locFile, result.Locus, result.Start, result.End, result.AltSeq)
						}

						// fmt.Printf(">%v %v(%v:%v) %v\n%v\n", locFile, result.locus, result.start, result.end, result.prod, result.altSeq)

					}
				} else {

					for _, file := range bsatstruct.ListOfFiles {
						if bsatseq.CheckVCFFormat(file) == true {
							result := bsatseq.Pos2Seq(*statInFile, file)
							locFile := strings.Replace(result.VcfFile, ".vcf", "", -1)
							fmt.Printf(">%v_%v(%v_%v)\n%v\n", locFile, result.Locus, result.Start, result.End, result.AltSeq)
						}
						// fmt.Printf(">%v %v(%v:%v) %v\n%v\n", locFile, result.locus, result.start, result.end, result.prod, result.altSeq)
					}
				}
			}

			// проверяет позицию из файла в снипах
		case "check_pos":

			if *statInFile != "" {
				// res := coord2gene(*statInFile)
				if *statVCF == "list" {

					/*
							go run bsatool.go  stat -b test_core -a check_pos -i uniq  --vcf=list --mkseq
							 go run bsatool.go  stat -b test_core -a check_pos -i uniq  --vcf=list --mkseq --reverse

						4013

					*/
					bsatf.CheckPosListVCF(*statInFile, *statMakeSeq)
				} else {

					bsatf.CheckPosList(*statInFile)
				}
			}

		case "Snp2SeqByLoc":
			/*
				go run bsatool.go stat  -b test_core -a Snp2SeqByLoc --vcf=list --locus=Rv1319c
			*/

			if *statLocus != "" && *statVCF == "list" {

				switch *statCircosTypeOf {
				case "":
					if len(*statInFile) != 0 {
						flist := bsatf.LoadFileNames(*statInFile)
						// fmt.Println(flist)
						for _, file := range flist {
							result := bsatseq.Snp2SeqByLoc(*statLocus, file)
							locFile := strings.Replace(result.VcfFile, ".vcf", "", -1)
							fmt.Printf(">%v\n%v\n", locFile, result.AltSeq)
							// fmt.Printf(">%v %v(%v:%v) %v\n%v\n", locFile, result.locus, result.start, result.end, result.prod, result.altSeq)

						}
					} else {
						for _, file := range bsatstruct.ListOfFiles {
							result := bsatseq.Snp2SeqByLoc(*statLocus, file)
							locFile := strings.Replace(result.VcfFile, ".vcf", "", -1)
							fmt.Printf(">%v\n%v\n", locFile, result.AltSeq)
							// fmt.Printf(">%v_%v(%v_%v)\n%v\n", locFile, result.locus, result.start, result.end, result.altSeq)
							// fmt.Printf(">%v %v(%v:%v) %v\n%v\n", locFile, result.locus, result.start, result.end, result.prod, result.altSeq)

						}
					}

				case "nc":
					bsatcalc.Locus2Matrix(*statLocus, bsatstruct.ListOfFiles, "nc")
				case "nc_coded":

					bsatcalc.Locus2Matrix(*statLocus, bsatstruct.ListOfFiles, "nc_coded")
				case "binary":

					bsatcalc.Locus2Matrix(*statLocus, bsatstruct.ListOfFiles, "binary")
				}

			}
		case "countSNP":
			/*
				go run bsatool.go stat  -b test_core -a countSNP --vcf=list
			*/

			if *statVCF == "list" {

				for _, file := range bsatstruct.ListOfFiles {
					fmt.Printf("#---- %v ------#\n", file)
					bsatvcf.CountSNP(file)

				}
			}
		}

	case "dev":

		bsatdb.LoadDB(*devDB)

		switch *devTask {
		case "dump":
			if *devPwd == "vsink" {
				fmt.Println(bsatstruct.FullGenesInfo)
				fmt.Println(bsatstruct.GenePositions)
			}
		}
		fmt.Println("DevMode")

	case "beta":

		bsatdb.LoadDB(*betaDB)

		switch *betaTask {

		case "complex":

			if *betaInFile != "" {
				bsatf.CalcCIdxFromFile(*betaInFile)
			}
		case "annotateFromList":
			if len(*betaInFile) != 0 {
				bsatf.AnnotFromFile(*betaInFile, bsatstruct.FullGenesInfo)
			}
		case "change":

			// if *betaInFile != "" {
			// bsatseq.MakeChangedSeq(2289241, 2288681, 2288681, "a", 10, 10)
			// }

		}

	case "parse_pileup":

		bsatdb.LoadDB(*annDB)

		if *pileupInFile != "" {
			bsatseq.Pileup2multifasta(*pileupInFile, *pileupTH, *pileupGName, *pileupGStrain, *pileupMkSeq, *pileupMinLen)
		}

	case "filter":
		bsatstruct.Flag.GbDP = *gbDP
		if *gbDP != 0 {

			for _, file := range bsatstruct.ListOfFiles {
				bsatf.FilterVcf(file, *gbDP)
			}
		}
	case "dpmap":
		bsatdb.LoadDB(*dpMapDB)
		allGenes := bsatservice.GetListOfGenes()
		for _, file := range bsatstruct.ListOfFiles {
			bsatf.DPmap(file, allGenes)
		}

	case "info":
		/*
		 go run bsatool.go info --db test_core --locus=Rv0019c --showas=gene
		 go run bsatool.go info --db test_core --locus=Rv0019c --showas=direct
		 go run bsatool.go info --db test_core --locus=Rv0278c --codons=1:30
		 go run bsatool.go info --db test_core --locus=Rv0278c --range=1:30
		 go run bsatool.go info --db test_core --locus=Rv2629 --showas=fasta --mask=191:a:c
		 go run bsatool.go info --db test_core --locus=Rv3915 --showas=fasta --mask=1040:a:c  --iupac
		 go run bsatool.go info --db test_core --locus=Rv0294 --showas=fasta --mask=726:t:c  --iupac --flank_right=400 --flank_left=200

		*/

		bsatdb.LoadDB(*infoDB)

		bsatseq.GetInfoByLocus()

	}

}
