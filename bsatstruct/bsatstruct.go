package bsatstruct

// import "bsatoolMod/src/bsatv

const (
	NCFlg = "NC"
	AAFlg = "AA"
)

var (
	GenomeSequence    []string // информация об генах, загруженная из базы
	ListOfFiles       []string
	FullGenesInfo     []GeneInfo
	GenePositions     = make(map[string]GCoords)
	GenomeDescription GenomeInfo
	TLMN              = "LMN"     // locus:Mutation:NAME
	TPMLN             = "PMLN"    // position:Mutation:locus:NAME
	TPMN              = "PMN"     // position:Mutation:NAME
	TLSAAN            = "LSAAN"   // locus:shortAA:codon:shortAA:name
	TLLAAN            = "LLAAN"   // locus:longAA:codon:longAA:name
	TLCN              = "LCN"     // locus:codon:name
	TSEN              = "SEN"     // start|end:name (any)
	TPNA              = "PNA"     // position~name~antibiotic
	TGLSAAAB          = "GLSAAAB" // geneOrLocus_shortAACodonshortAA whiB6_P38L_SM  Rv1258c_V219A_SM

	/*

		LmnRegExp  = `^(\w+)\W+(\d+)(\D)>(\D)\s+(.*)` // LMN Regular expression
			PmlnRegExp = `^(\d+)_(\D)>(\D)\W+(\w+)\W+(.*)`
			PmnRegExp  = `(\d+)_(\D)>(\D)\W+(.*)$`

			LsaanRegExp    = `^(\w+)\W+(\D{1})(\d+)(\D{1})\W+(.*)`
			LlaanRegExp    = `^(\w+)\W+(\w{3})\W+(\d+)\W+(\w{3})\W+(.*)`
			LcnRegExp      = `^(\w+)\W+codon(\d+)\W+(.*)`
			SenRegExp      = `^(\d+)\|(\d+):(\w+)`
			TpnaRegExp     = `^(\d+)~(\S+)~(\w+)`
			GlsaaABRegExp  = `^(\S+)_(\w)(\d+)(\w)_(\w+)`
			ChkPosRegExp   = `^check\W+(\d+)\W+(.*)`
			ChkCodonRegExp = `^check\W+(\w+)\W+(\d+)\W+(.*)`

	*/

	SnpCacheMap = make(map[string][]SnpInfo)
	Flag        Flags
	RulesArr    []RulesInfo
	ResStart    int
	ResEnd      int
)

type (
	// GeneInfo 'GeneInfo ....'

	GeneInfo struct {
		Locus, Name, Product, Direction, GeneID, ProteinID, Note, GOA, TypeOf string
		Start, End, NucCore                                                   int
		PDB, InterPro, ProSite                                                []string
	}
	GenomeInfo struct {
		Organism, Strain, Version string
		Start, End, NucCore       int
	}
	GCoords struct {
		Start, End int    // начало и конец гена
		Type       string // тип гена CDS
	}

	SnpCheckInfo struct {
		Locus, PosInGene, CodonNbrInG, Ref, Alt, Name, TypeOf, APos, AASref, AASalt, AALref, AALalt, Raw, Tag, AB string
		StartRange, EndRange                                                                                      int
	}

	SnpInfo struct {
		/*
					APos абсолютная позиция в геноме
						PosInGene позиция в гене
			PosInCodonG позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)

		*/
		APos, PosInGene, PosInCodonG, CodonNbrInG, GeneLen, Start, End, NucCore, TangIdxVal, DP, Indel int
		RefCodon, AltCodon, RefAA, AltAA, Locus,
		Direction, NucInPosCoding, NucInGenomeRef, NucInGenomeAlt, Product, Name,
		RefAAShort, AltAAShort, Mutation, Tang, Alt, Note, ReportType, ProteinID, GeneID, GOA, TiTv, TypeOf, ComplexIndex, FName, IndelType, IndelRef,
		IndelAlt string
		InterPro, PDB, ProSite []string
	}
	CheckSNP struct {
		FileName, FoundSNP, ExactWithRule, RuleNames string
	}
	RulesInfo struct {
		Name     string
		Variants []string
		Lenght   int
	}

	AllPositionsInGene struct {
		Pos             int
		Alt, Ref, Locus string
	}

	GC3Type struct {
		GC3Alt, GC3Ref, GCalt, GCref float64
		Locus                        string
	}

	// JaroWinklerInfo is....
	JaroWinklerInfo struct {
		Locus, Product string
		JWDist         float64
	}

	StatInfo struct {
		Pos, Count, Perc        int
		FilesWith, FilesWithout string
	}

	SeqInfo struct {
		Name, Seq, TypeOfSeq string //
		UsedPositions        []int
		DebugInfo            []string
		// Len       int
	}

	Flags struct {
		// StatTH, StatDP, GbDP, StatFlankLeft, StatFlankRight                                                int
		// AnnShowFileName, GbDebug, GbIndex, GbVerbose, GbInDel, GbIGR, StatAsGene, StatAsTable, StatReverse,AnnMakeSeqRef,BbRandomize bool
		// AnnBench,StatInRule, GbPort string

		GbWeb     bool
		GbXLSX    string
		GbVerbose bool
		GbIndex   bool
		GbPort    string
		GbNoSeq   bool
		GbDebug   bool
		// gbLog            = kingpin.Flag("log", "write log file").Default("false").Bool()
		GbExcludeGenes   string
		GbExcludeRegions string
		GbExcludeSnp     string
		GbRandomize      bool
		DbName           string
		DbGenbank        string
		GbDP             int
		GbIGR            bool
		GbInDel          bool

		AnnDB            string
		AnnVCF           string
		AnnMakeSeq       string
		AnnOutFormat     string
		AnnMakeSeqRef    bool
		AnnWithFilenames bool
		AnnBench         string
		AnnSeqLen        int
		AnnShowFileName  bool
		AnnPosFile       string
		AnnMinPosNbr     int
		AnnAltRange      string

		StatDB         string
		StatTask       string
		StatInFile     string
		StatInFileList string
		StatOutFile    string
		StatTypeOf     string
		StatInRule     string
		StatAll        bool
		StatBench      string
		StatMakeSeq    bool
		// statRefGenomeName   = statAction.Flag("ref", "set custom genome name").String()
		StatCircosTypeOf    string
		StatBedTypeOf       string
		StatCircosBandColor string
		StatGenomeName      string
		StatShowAnnotation  bool
		StatNbrOfSNP        int
		StatMinLocusCount   int
		StatGroupFromFile   string
		StatVCF             string
		StatLocus           string
		StatDnDsCountTh     int
		StatMkSeq           bool
		StatCircos          bool
		StatCircosGC        bool
		StatFlankLeft       int
		StatFlankRight      int
		StatAsTable         bool
		StatAsGene          bool
		StatReverse         bool
		StatTH              int
		StatDnDsStat        bool
		InfoLocus           string
		InfoDB              string
		InfoRanges          string
		InfoCodons          string
		InfoShowAs          string
		InfoMask            string
		InfoFlankLeft       int
		InfoFlankRight      int
		InfoIUPAc           bool
		DevDB               string
		DevTask             string
		DevPwd              string
		BetaTask            string
		BetaDB              string
		BetaInFile          string

		PileupInFile  string
		PileupTH      int
		PileupGStrain string
		PileupMkSeq   bool
		PileupShowSeq bool
		PileupGName   string
		PileupCircos  bool
		PileupDB      string
		PileupMinLen  int
	}
	AltStringResult struct {
		Start, End                   int
		Locus, Prod, AltSeq, VcfFile string
	}
	RangePosInfo struct {
		Start, End, Len, Doubles          int
		Seq, Prod, GeneName, Note, Genome string
		GC                                float64
	}

	RangeArray struct {
		Start, End int
	}
	GenomeMapInfo struct {
		Start, End             int
		Locus, Product, TypeOf string
	}
)
