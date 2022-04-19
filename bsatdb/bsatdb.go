package bsatdb

/*

Модуль для работы с базой данных. Ее создание. Загрузка.


*/

import (
	"bsatoolMod/src/bsatstruct"
	"bufio"
	"compress/lzw"
	"encoding/gob"
	"fmt"
	"io"
	"log"
	// "bytes"
	// "fmt"
	// "strings"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

var (
	logoLoc = `
 _______    ________     __  ___________  ______      ______    ___       
|   _  "\  /"       )   /""\("     _   ")/    " \    /    " \  |"  |      
(. |_)  :)(:   \___/   /    \)__/  \\__/// ____  \  // ____  \ ||  |      
|:     \/  \___  \    /' /\  \  \\_ /  /  /    ) :)/  /    ) :)|:  |      
(|  _  \\   __/  \\  //  __'  \ |.  | (: (____/ //(: (____/ //  \  |___   
|: |_)  :) /" \   :)/   /  \\  \\:  |  \        /  \        /  ( \_|:  \  
(_______/ (_______/(___/    \___)\__|   \"_____/    \"_____/    \_______) 
` +
		`BSATool - Bacterial Snp Annotation Tool ` + "\n" +
		`      Laboratory of Social and Epidemic Infections
 Scientific Centre for Family Health and Human Reproduction Problems
     	(c) V.Sinkov, Irkutsk, Russia, 2017-2022                                   
                                                  
	`
	genomeSeqSlice  []string // информация об генах, загруженная из базы
	gInfo           bsatstruct.GenomeInfo
	geneCoordinates = make(map[string]bsatstruct.GCoords)
	allGenesVal     []bsatstruct.GeneInfo
	// GeneIn
)

func readDB(file string) []bsatstruct.GeneInfo {
	var (
		gene []bsatstruct.GeneInfo
	)

	gobFile, err := os.Open(file)
	if err != nil {
		log.Println(file, err.Error())
	}

	defer func(gobFile *os.File) {
		err := gobFile.Close()
		if err != nil {
			log.Println(file, err.Error())
		}
	}(gobFile)

	compressedGobFile := lzw.NewReader(gobFile, lzw.LSB, 8)
	defer func(compressedGobFile io.ReadCloser) {
		err := compressedGobFile.Close()
		if err != nil {
			log.Println(file, err.Error())
		}
	}(compressedGobFile)
	gobParser := gob.NewDecoder(compressedGobFile)
	_ = gobParser.Decode(&gene)
	_ = gobParser.Decode(&genomeSeqSlice)
	_ = gobParser.Decode(&gInfo)
	_ = gobParser.Decode(&geneCoordinates)

	return gene

}

func LoadDB(file string) {
	locGene := readDB(file)
	bsatstruct.FullGenesInfo = locGene
	bsatstruct.GenomeSequence = genomeSeqSlice
	bsatstruct.GenePositions = geneCoordinates
	bsatstruct.GenomeDescription = gInfo

}

func parseGeneBankFile(file string) (GenomeData []bsatstruct.GeneInfo, GenomeSplice []string) {
	// функция для считывания файла в формате генбанка и занесения строк в массив linesFromGB

	var (
		qenomeSeq    = regexp.MustCompile(`\s+(\d+)\s+(\w{6}.*)`)
		regDelSpaces = regexp.MustCompile(`(\s+)(.*)`)
		checkOrigin  = regexp.MustCompile(`(\d+\s\w{10}\s\w{10}\s\w{10})`)

		startOrigin    = regexp.MustCompile(`(ORIGIN\s+)|(translation=)`)
		feature        = regexp.MustCompile(`^\s{5}(\w+)\s+complement\W(\d)+\W+(\d)+\W|^\s{5}(\w+)\s+(\d)+\W+(\d+)`)
		makeAnchors    = regexp.MustCompile(`^(/)`)
		genomeSource   = regexp.MustCompile(`organism=\W(.*)\W`)
		genomeStrain   = regexp.MustCompile(`strain=\W(.*)\W`)
		genomeVersion  = regexp.MustCompile(`^VERSION\s+(.*)`)
		genomeStartEnd = regexp.MustCompile(`source\W+(\d+)..(\d+)`)

		gStart, gEnd int
		gVersion     string

		resString []string

		extractedData, extractedDataFeat []string

		originBlock, CDSblock, firstCDS int

		changedStr, organismName, organismStrain string
		noteBuffer                               strings.Builder
	)

	fmt.Println(logoLoc)
	f, err := os.Open(file) // открываем файл

	if err != nil {
		fmt.Println(err)

	}
	defer func(f *os.File) {
		err := f.Close()
		if err != nil {
			log.Println(file, err.Error())
		}
	}(f)
	extractedData = append(extractedData, "START_BLOCK") // начало блока записи
	extractedDataFeat = append(extractedDataFeat, "START_BLOCK")
	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()
		// размер генома
		for _, gStartEnd := range genomeStartEnd.FindAllStringSubmatch(scanTxt, -1) {
			gStart, _ = strconv.Atoi(gStartEnd[1])
			gEnd, _ = strconv.Atoi(gStartEnd[2])
		}
		//  организм
		for _, gName := range genomeSource.FindAllStringSubmatch(scanTxt, -1) {
			organismName = gName[1]
			fmt.Printf("Organism:%v\n", organismName)
		}

		// штамм
		for _, gStrain := range genomeStrain.FindAllStringSubmatch(scanTxt, -1) {
			organismStrain = gStrain[1]
			fmt.Printf("Strain:%v\n%v%v...\r", organismStrain, "Working on ", file)
		}
		for _, genomevesion := range genomeVersion.FindAllStringSubmatch(scanTxt, -1) {
			gVersion = genomevesion[1]
			fmt.Printf("Version:%v\n", gVersion)
		}

		for _, sOrigin := range startOrigin.FindAllStringSubmatch(scanTxt, -1) {
			if len(sOrigin[1]) != 0 || len(sOrigin[2]) != 0 {
				CDSblock = 0
			}
		}

		for _, rfeat := range feature.FindAllStringSubmatch(scanTxt, -1) {
			if rfeat[1] != "gene" && rfeat[4] != "gene" && strings.Contains(rfeat[0], "source") == false {
				CDSblock = 1

				if rfeat[1] != "" {
					scanTxt = strings.Replace(scanTxt, rfeat[1], fmt.Sprintf("*%v", rfeat[1]), -1)
				} else if rfeat[4] != "" {
					scanTxt = strings.Replace(scanTxt, rfeat[4], fmt.Sprintf("*%v", rfeat[4]), -1)

				}

			} else if rfeat[1] == "gene" || rfeat[4] == "gene" && strings.Contains(rfeat[0], "source") == false {

				CDSblock = 0

			}
		}

		for _, origninFound := range checkOrigin.FindAllStringSubmatch(scanTxt, -1) {
			if len(origninFound[1]) != 0 {
				originBlock = 1
				CDSblock = 0

			} else {
				originBlock = 0
			}
		}

		if CDSblock == 1 {
			if firstCDS == 0 {

				firstCDS = 1

			}
			changedStr = regDelSpaces.ReplaceAllString(scanTxt, "$2 ") // удаляем пробелы

			changedStr = strings.Replace(changedStr, "/note=", "!note=", -1) // меняем / на ! для дальнейшего парсинга

			changedStr = strings.Replace(changedStr, "/product=", "!product=", -1) // см выше.

			changedStr = strings.Replace(changedStr, "\"", "", -1)

			changedStr = makeAnchors.ReplaceAllString(changedStr, "!!")
			if strings.Index(changedStr, "!!") != 0 && strings.Index(changedStr, "*") != 0 {

				noteBuffer.WriteString(changedStr)

			} else {
				if noteBuffer.Len() != 0 {

					extractedData = append(extractedData, strings.TrimSpace(strings.TrimPrefix(strings.Replace(noteBuffer.String(), "!", "\n!", -1), "!")))
					noteBuffer.Reset()

				}

			}

			if strings.Index(changedStr, "!!") == 0 {

				extractedData = append(extractedData, strings.TrimSpace(strings.TrimPrefix(changedStr, "!!")))

			} else if strings.Index(changedStr, "*") == 0 {

				if firstCDS == 1 {

					extractedData = append(extractedData, "START_OF")
					firstCDS = 2
				} else if firstCDS == 2 {
					extractedData = append(extractedData, "END_OF")
					extractedData = append(extractedData, "START_OF")
				}
				extractedData = append(extractedData, strings.TrimPrefix(changedStr, "*"))

			}

		}

		if originBlock == 1 {

			for _, genomeMatch := range qenomeSeq.FindAllStringSubmatch(scanTxt, -1) {

				resString = append(resString, strings.Replace(genomeMatch[2], " ", "", -1))
			}
		}

	}
	extractedData = append(extractedData, "END_OF")
	extractedData = append(extractedData, "END_BLOCK")
	GenomeSplice = strings.SplitAfter(strings.Join(resString, ""), "")

	var (
		lDir, lLoc, lName, lProd, lNote, gID, pID, lGOA, lType string
		lPDBArr, lInterProArr, lProSiteArr                     []string
		nucCore, igrCount                                      int
		cdsStEnd                                               = regexp.MustCompile(`^(\w+)\s+(\d+)\W+(\d+)|^(\w+)\s+complement\W(\d+)\W+(\d+)`)
		cdsLocus                                               = regexp.MustCompile(`locus_tag=(.*)`)
		cdsName                                                = regexp.MustCompile(`gene=(.*)`)
		cdsProd                                                = regexp.MustCompile(`product=(.*)`)
		cdsNote                                                = regexp.MustCompile(`note=(.*)`)
		cdsgID                                                 = regexp.MustCompile(`db_xref=GeneID:(\d+)`)
		cdsprotID                                              = regexp.MustCompile(`protein_id=(.*)`)
		cdsGOA                                                 = regexp.MustCompile(`db_xref=GOA:(.*)`)
		cdsPDB                                                 = regexp.MustCompile(`db_xref=PDB:(.*)`)
		cdsInterPro                                            = regexp.MustCompile(`db_xref=InterPro:(.*)`)
		cdsProSite                                             = regexp.MustCompile(`inference=protein motif:PROSITE:(.*)`)

		startOfBlock, endOfBlock, cdsStart, cdsEnd, lStart, lEnd int

		cdsCount       = map[string]int{}
		igensS, igensE []int

		leftGene []string
	)

	for _, val := range extractedData {
		// fmt.Println(val)
		if strings.Index(val, "START_BLOCK") != -1 {
			startOfBlock = 1
			endOfBlock = 0
		} else if strings.Index(val, "END_BLOCK") != -1 {
			endOfBlock = 1
			startOfBlock = 0
		}
		if strings.Index(val, "START_OF") != -1 {
			cdsStart = 1
			cdsEnd = 0
			// cdsCount++
		} else if strings.Index(val, "END_OF") != -1 {
			cdsStart = 0
			cdsEnd = 1
		}

		if startOfBlock == 1 && endOfBlock == 0 {

			if cdsStart == 1 && cdsEnd == 0 {

				for _, cdsStEndMatch := range cdsStEnd.FindAllStringSubmatch(val, -1) {

					if cdsStEndMatch[2] != "" {

						lType = cdsStEndMatch[1]
						lStart, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[2]))
						lEnd, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[3]))

						lDir = "f"

						cdsCount[cdsStEndMatch[1]] = cdsCount[cdsStEndMatch[1]] + 1

					} else if cdsStEndMatch[5] != "" {

						lType = cdsStEndMatch[4]
						lStart, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[5]))
						lEnd, _ = strconv.Atoi(strings.TrimSpace(cdsStEndMatch[6]))

						lDir = "r"

						cdsCount[cdsStEndMatch[4]] = cdsCount[cdsStEndMatch[4]] + 1

					}

				}

				for _, cdsNameMatch := range cdsName.FindAllStringSubmatch(val, -1) {
					lName = strings.TrimSpace(strings.Replace(cdsNameMatch[1], " ", "", -1))

				}

				for _, cdsLocusMatch := range cdsLocus.FindAllStringSubmatch(val, -1) {
					lLoc = strings.TrimSpace(strings.Replace(cdsLocusMatch[1], " ", "", -1))

				}
				for _, cdsProdMatch := range cdsProd.FindAllStringSubmatch(val, -1) {
					lProd = strings.TrimSpace(cdsProdMatch[1])

				}
				for _, cdsNoteMatch := range cdsNote.FindAllStringSubmatch(val, -1) {
					lNote = strings.TrimSpace(cdsNoteMatch[1])

				}
				for _, cdsgIDMatch := range cdsgID.FindAllStringSubmatch(val, -1) {
					gID = strings.TrimSpace(strings.Replace(cdsgIDMatch[1], " ", "", -1))

				}
				for _, cdsgprotIDMatch := range cdsprotID.FindAllStringSubmatch(val, -1) {
					pID = strings.TrimSpace(strings.Replace(cdsgprotIDMatch[1], " ", "", -1))

				}
				for _, cdsGOAMatch := range cdsGOA.FindAllStringSubmatch(val, -1) {
					lGOA = strings.TrimSpace(strings.Replace(cdsGOAMatch[1], " ", "", -1))

				}
				for _, cdsPDBMatch := range cdsPDB.FindAllStringSubmatch(val, -1) {
					lPDBArr = append(lPDBArr, strings.TrimSpace(strings.Replace(cdsPDBMatch[1], " ", "", -1)))

				}
				for _, cdsInterProMatch := range cdsInterPro.FindAllStringSubmatch(val, -1) {
					lInterProArr = append(lInterProArr, strings.TrimSpace(strings.Replace(cdsInterProMatch[1], " ", "", -1)))

				}
				for _, cdsProSiteMatch := range cdsProSite.FindAllStringSubmatch(val, -1) {
					lProSiteArr = append(lProSiteArr, strings.TrimSpace(strings.Replace(cdsProSiteMatch[1], " ", "", -1)))

				}

				if gID != "" && lGOA == "" {
					nucCore = 1

				} else if gID == "" && lGOA != "" {
					nucCore = 1

				} else if gID == "" && lGOA == "" {
					nucCore = 2
				}

			} else if cdsStart == 0 && cdsEnd == 1 {

				if lLoc == "" && lName == "" {
					lLoc = fmt.Sprintf("%v_%v_%v", lType, lStart, lEnd)

				} else if lLoc == "" && lName != "" {
					lLoc = lName
				}

				if lProd == "" && lNote != "" {
					lProd = lNote
				} else if lProd == "" && lNote == "" {
					lProd = lType
				}

				geneInfo := bsatstruct.GeneInfo{Start: lStart, End: lEnd, Locus: lLoc, Direction: lDir,
					Product: lProd, Name: lName, GeneID: gID, ProteinID: pID, Note: lNote, GOA: lGOA, TypeOf: lType, PDB: lPDBArr, InterPro: lInterProArr, ProSite: lProSiteArr}

				GenomeData = append(GenomeData, geneInfo)

				if lType == "CDS" {

					geneCoordinates[lLoc] = bsatstruct.GCoords{Start: lStart, End: lEnd, Type: "CDS"}

				}

				lName, lStart, lEnd, lLoc, lDir, lProd, gID, pID, lGOA, lNote, lType = "", 0, 0, "", "", "", "", "", "", "", ""
				lPDBArr = nil
				lInterProArr = nil
				lProSiteArr = nil

			}
		}

	}

	for key, val := range cdsCount {
		fmt.Printf("\nFound: %v %v", val, key)
	}
	if nucCore == 0 {
		fmt.Println("\nGene Info: Genbank Id Nomenclature", nucCore)
	} else if nucCore == 1 {
		fmt.Println("\nGene Info: high-quality Gene Ontology (GO) annotations", nucCore)
	} else if nucCore == 2 {
		fmt.Println("\nGene Info: unknown", nucCore)
	}

	gInfo = bsatstruct.GenomeInfo{NucCore: nucCore, Organism: organismName, Start: gStart, End: gEnd, Strain: organismStrain, Version: gVersion}

	sort.Slice(
		GenomeData, func(i, j int) bool {
			return GenomeData[i].Start < GenomeData[j].Start
		})

	for i, gene := range GenomeData {

		igensS = append(igensS, gene.Start)
		igensE = append(igensE, gene.End)
		leftGene = append(leftGene, gene.Locus)

		if i >= 1 && i < len(GenomeData)-1 {

			igrCount++
			checkEndOfCDS := (igensE[i-1] + 1) - gene.Start

			if checkEndOfCDS < 0 {
				igr := bsatstruct.GeneInfo{Start: igensE[i-1] + 1, End: igensS[i] - 1, Locus: fmt.Sprintf(
					"IGR_%v_%v", igensE[i-1]+1, igensS[i]-1), Direction: "f",
					Product: "Intergenic region (IGR)", Name: "IGR", TypeOf: "IGR"}
				GenomeData = append(GenomeData, igr)

				geneCoordinates[fmt.Sprintf(
					"IGR_%v_%v", igensE[i-1]+1, igensS[i]-1)] = bsatstruct.GCoords{Start: igensE[i-1] + 1, End: igensS[i] - 1, Type: "IGR"}
			}

		}

	}
	fmt.Println("Found IGR regions:", igrCount)

	sort.Slice(
		GenomeData, func(i, j int) bool {
			return GenomeData[i].Start < GenomeData[j].Start
		})

	return GenomeData, GenomeSplice

}

func writeDB(file string, gene []bsatstruct.GeneInfo) {

	gobFile, err := os.Create(file)
	if err != nil {
		log.Println(err.Error())
	}
	if gobFile != nil {
		defer func(gobFile *os.File) {
			err := gobFile.Close()
			if err != nil {
				log.Println(file, err.Error())
			}
		}(gobFile)
	}

	compressedGobFile := lzw.NewWriter(gobFile, lzw.LSB, 8)
	defer func(compressedGobFile io.WriteCloser) {
		err := compressedGobFile.Close()
		if err != nil {
			log.Println(file, err.Error())
		}
	}(compressedGobFile)
	gobParser := gob.NewEncoder(compressedGobFile)
	_ = gobParser.Encode(&gene)
	_ = gobParser.Encode(&genomeSeqSlice)
	_ = gobParser.Encode(&gInfo)
	_ = gobParser.Encode(&geneCoordinates)

	fmt.Println(file, " was created successfully.")

}

func CreateDB(gbFile string, dbName string) {
	allGenesVal, genomeSeqSlice = parseGeneBankFile(gbFile)

	writeDB(dbName, allGenesVal)
}
