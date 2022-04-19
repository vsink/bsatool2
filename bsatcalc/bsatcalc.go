package bsatcalc

import (
	"bsatoolMod/src/amino"
	"bsatoolMod/src/bsatf"
	"bsatoolMod/src/bsatseq"
	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bsatoolMod/src/bsatvcf"
	codon "bsatoolMod/src/codonPkg"
	"bufio"
	"fmt"
	"github.com/xrash/smetrics"
	"log"
	"math/rand"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"
)

func CalcDnDs(vcfFile string) {

	var (
		allLocuses []string

		snps         []bsatstruct.SnpInfo
		altPositions = make(map[string][]bsatstruct.AllPositionsInGene)
	)

	for key, val := range bsatstruct.GenePositions {
		if val.Type == "CDS" {
			allLocuses = append(allLocuses, key)
		}
	}

	sort.Strings(allLocuses)

	snpsChan := make(chan []bsatstruct.SnpInfo)

	go func() {
		snpsChan <- bsatvcf.MakeSnps(vcfFile)
	}()
	snps = <-snpsChan

	for _, val := range snps {
		if val.TypeOf == "CDS" && bsatseq.ChkPosExists(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

			altPositions[val.Locus] = append(
				altPositions[val.Locus], bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})
		}

	}

	for _, allloc := range allLocuses {

		if len(altPositions[allloc]) > 2 {
			dndsChan := make(chan []string)
			go func() {

				dndsChan <- GetDnDsByLocus(allloc, altPositions[allloc])
			}()
			dndsRes, ok := <-dndsChan

			if ok {

				fmt.Println(fmt.Sprintf("%v(%v)\t", allloc, bsatseq.GetProdByName(allloc)), dndsRes[1])
				close(dndsChan)
			}

		} else {

		}

	}

}

func CalcGC3Val(snps []bsatstruct.SnpInfo) []bsatstruct.GC3Type {

	var (
		altPositions = make(map[string][]bsatstruct.AllPositionsInGene)

		gcArray []bsatstruct.GC3Type

		allLocusUnsort, allLocuses []string
	)

	for _, val := range snps {

		if val.TypeOf == "CDS" {

			altPositions[val.Locus] = append(
				altPositions[val.Locus], bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})

			allLocusUnsort = append(allLocusUnsort, val.Locus)
		}

	}

	allLocuses = bsatservice.RmStrDoubles(allLocusUnsort)

	for _, loc := range allLocuses {

		refS := bsatseq.GetGeneSeq(loc)

		altS := bsatseq.MakeAltString(loc, altPositions[loc])

		gcAlt, _, _, gc3Alt := codon.GcCodonCalc(altS)
		gcRef, _, _, gc3Ref := codon.GcCodonCalc(refS)

		gcArray = append(gcArray, bsatstruct.GC3Type{Locus: loc, GC3Alt: gc3Alt, GC3Ref: gc3Ref, GCalt: gcAlt, GCref: gcRef})

	}

	return gcArray

}

func CalcJaroWinklerDist(file string, print bool) []bsatstruct.JaroWinklerInfo {

	var (
		altPositions = make(map[string][]bsatstruct.AllPositionsInGene)
		validData    []string
		altSequences []string
		jwRes        []bsatstruct.JaroWinklerInfo

		jarwinkl float64
	)

	if len(bsatstruct.SnpCacheMap) == 0 {
		qSNP := &bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery)}
		go qSNP.Request()

		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
	}

	for fname, snps := range bsatstruct.SnpCacheMap {
		if fname == file {
			for _, val := range snps {

				if val.TypeOf == "CDS" {
					altPositions[val.Locus] = append(
						altPositions[val.Locus], bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})
				}
			}
		}
	}

	// выбор генов, в которых обнаружены мутации
	for key, val := range altPositions {

		if len(val) > 1 {
			validData = append(validData, key)
		}
	}

	sort.Strings(validData)

	for _, val := range validData {

		altS := bsatseq.MakeAltString(val, altPositions[val])
		altSequences = append(altSequences, altS)
		refS := bsatseq.GetGeneSeq(val)

		jarwinkl = smetrics.JaroWinkler(refS, altS, 0.7, 4)
		jwRes = append(jwRes, bsatstruct.JaroWinklerInfo{Locus: val, JWDist: jarwinkl})

		if print == true {
			fmt.Printf("L:%v JW:%v\n", val, jarwinkl)
		}

	}

	return jwRes

}

func MakeSNPStatistic() {
	// var countSNPs = 1
	var (
		pos                  = make(map[int]int)
		alt                  = make(map[int]bsatstruct.SnpInfo)
		f                    = make(map[int][]string)
		positions, posUnsort []int
	)

	files := &bsatstruct.ListOfFiles

	upperLimit := len(*files)

	for i, file := range *files {

		fmt.Printf("Reading: %v (%v from %v)%v \r", file, i+1, len(*files), strings.Repeat(" ", 60))

		qSNP := &bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery)}
		go qSNP.Request()

		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
	}

	for fname, snps := range bsatstruct.SnpCacheMap {

		for _, val := range snps {
			if pos[val.APos] <= upperLimit {
				pos[val.APos] = pos[val.APos] + 1        // count
				alt[val.APos] = val                      // pos
				f[val.APos] = append(f[val.APos], fname) // files
				posUnsort = append(posUnsort, val.APos)  // array of positions

			}

		}
		positions = bsatservice.GenUnique(posUnsort)
		sort.Ints(positions)

	}
	var stat []bsatstruct.StatInfo

	for _, p := range positions {
		perc := (pos[p] * 100) / upperLimit
		filesNotInList := bsatservice.CompareSlices(*files, f[p])
		stat = append(
			stat, bsatstruct.StatInfo{Pos: p, Count: pos[p], Perc: perc, FilesWith: strings.Join(f[p], ",\n"), FilesWithout: strings.Join(
				filesNotInList, ","+
					"\n")})

	}
	bsatservice.PrintStatsInWeb(stat, bsatstruct.Flag.GbPort)

}

// func MakeMatrix(typeof string, fileOut string, verbose bool) {
//
// 	var (
// 		AllPos     []int
// 		allLocuses []string
// 		buffer     strings.Builder
// 		headers    strings.Builder
// 		posCount   = make(map[int]int)
// 	)
//
// 	files := &bsatstruct.ListOfFiles
//
// 	for i, file := range *files {
// 		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
// 		qSNP := &bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery)}
// 		go qSNP.Request()
// 		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
// 		snpRes := <-qSNP.OutChan
// 		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
// 		if verbose == true {
// 			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
// 		}
//
// 	}
//
// 	for key, val := range geneCoordinates {
// 		if val.Type == "CDS" {
// 			allLocuses = append(allLocuses, key)
// 		}
// 	}
// 	sort.Strings(allLocuses)
//
// 	var (
// 		exGenes = make(map[int]int)
// 		exSNP   = make(map[int]int)
// 	)
// 	if *gbExcludeGenes != "" {
// 		exGenes = bsatf.LoadExclGenes(*gbExcludeGenes)
// 	} else if *gbExcludeSnp != "" {
// 		exSNP = bsatf.LoadExclSNP(*gbExcludeSnp)
// 	}
//
// 	switch typeof {
//
// 	case "binary":
//
// 		matrixBinary(fileOut, exGenes, exSNP)
//
// 	// case "table":
//
// 	// 	matrixTable(fileOut)
//
// 	case "nc":
// 		var posFreq = map[int][]string{}
// 		pos := make(map[int]string)
//
// 		headers.WriteString("Pos\tRef\t")
// 		// allPosChan := make(chan []int)
// 		// go func() {
// 		// 	allPosChan <- getAllPosFromCacheMap()
// 		// }()
// 		// AllPos = <-allPosChan
// 		// for i, file := range files {
// 		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
// 		for fname, snps := range bsatstruct.SnpCacheMap {
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))
//
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
//
// 			// snps := parserVCF(file, false, allGenesVal)
//
// 			for _, val := range snps {
//
// 				pos[val.APos] = val.Alt
// 				// AllPosUnsort = append(AllPosUnsort, val.APos)
//
// 			}
// 			for _, allpos := range AllPos {
//
// 				if pos[allpos] != "" {
// 					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", pos[allpos]))
// 					// fmt.Println(posFreq[allpos])
//
// 				} else {
// 					posFreq[allpos] = append(posFreq[allpos], bsatseq.GetNucFromGenomePos(allpos))
//
// 				}
//
// 			}
//
// 		}
// 		// AllPos = unique(AllPosUnsort)
// 		// // allLocuses = removeStringDuplicates(allLocusUnsort)
// 		// sort.Ints(AllPos)
//
// 		for _, allpos := range AllPos {
// 			// if buffer.Len() == 0 {
// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", bsatseq.GetNucFromGenomePos(allpos), "\t", strings.Join(posFreq[allpos], "\t")))
// 			// fmt.Println(buffer.String())
//
// 		}
// 		headers.WriteString("\n")
// 		matrixPrint(headers, buffer, fileOut)
// 	case "nc_coded":
// 		// A, T, G,C, ---> 0, 1, 2, 3
// 		var posFreq = map[int][]string{}
// 		pos := make(map[int]string)
// 		posCount := make(map[int]int)
// 		files := &bsatstruct.ListOfFiles
// 		headers.WriteString("Pos\t")
// 		// allPosChan := make(chan []int)
// 		// go func() {
// 		// 	allPosChan <- getAllPosFromCacheMap()
// 		// }()
// 		// AllPos = <-allPosChan
// 		// for i, file := range files {
// 		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
// 		for fname, snps := range bsatstruct.SnpCacheMap {
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))
//
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
//
// 			// snps := parserVCF(file, false, allGenesVal)
//
// 			for _, val := range snps {
//
// 				pos[val.APos] = val.Alt
// 				// AllPosUnsort = append(AllPosUnsort, val.APos)
// 				posCount[val.APos] = posCount[val.APos] + 1
// 			}
// 			for _, allpos := range AllPos {
//
// 				if pos[allpos] != "" {
//
// 					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", nuc2IntCode(pos[allpos])))
// 					// fmt.Println(posFreq[allpos])
//
// 				} else {
// 					posFreq[allpos] = append(posFreq[allpos], nuc2IntCode(bsatseq.GetNucFromGenomePos(allpos)))
//
// 				}
//
// 			}
//
// 		}
// 		// AllPos = unique(AllPosUnsort)
// 		// // allLocuses = removeStringDuplicates(allLocusUnsort)
// 		// sort.Ints(AllPos)
//
// 		for _, allpos := range AllPos {
// 			// if buffer.Len() == 0 {
// 			if posCount[allpos] <= len(*files)-1 {
// 				buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))
// 				// fmt.Println(buffer.String())
// 			}
//
// 		}
// 		headers.WriteString("\n")
// 		matrixPrint(headers, buffer, fileOut)
//
// 	case "summary":
// 		var (
// 			posFreq = map[int][]string{}
// 			// pos                    = make(map[int]string)
// 			groupRegexp            = regexp.MustCompile(`^(\S+)\W+(\w+)\W+(\w+)`)
// 			group                  = make(map[string]string)
// 			label                  = make(map[string]string)
// 			groupHeader, groupBody string
// 			fileGroup              = make(map[int][]string)
// 			fileLabel              = make(map[int][]string)
// 			PosInGenome            = make(map[string]map[int]string)
// 			PosAdditionalInfo      = make(map[int]string)
// 			genomes                []string
// 		)
//
// 		if *statGroupFromFile != "" {
// 			f, err := os.Open(*statGroupFromFile) // открываем файл
//
// 			if err != nil {
// 				fmt.Println(err)
//
// 			}
// 			defer f.Close()
//
// 			scanner := bufio.NewScanner(f) //  новый сканер
//
// 			for scanner.Scan() {
//
// 				scanTxt := scanner.Text()
//
// 				for _, grpVal := range groupRegexp.FindAllStringSubmatch(scanTxt, -1) {
// 					group[strings.ToUpper(grpVal[1])] = strings.ToUpper(grpVal[2])
//
// 					label[strings.ToUpper(grpVal[1])] = grpVal[3]
//
// 					// fmt.Println(grpVal)
// 				}
// 			}
//
// 		}
//
// 		// fmt.Println(group)
// 		if len(group) != 0 {
// 			groupHeader = "Group\t"
// 			if len(label) != 0 {
// 				groupHeader = fmt.Sprintf("%vLabel\t", groupHeader)
// 			}
// 		}
//
// 		headers.WriteString(fmt.Sprintf("Pos\tLocus\t%vGene\tAAchange\tDirection\tMutation\tCIndex\tCIndex_Res\tPosInCodon\t", groupHeader))
// 		// allPosChan := make(chan []int)
// 		// go func() {
// 		// 	allPosChan <- getAllPosFromCacheMap()
// 		// }()
// 		// AllPos = <-allPosChan
// 		// for i, file := range files {
// 		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
//
// 		for fname, snps := range bsatstruct.SnpCacheMap {
// 			PosInGenome[fname] = make(map[int]string)
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))
//
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
// 			genomes = append(genomes, fname)
// 			// snps := parserVCF(file, false, allGenesVal)
//
// 			for _, val := range snps {
//
// 				// pos[val.APos] = val.Alt
//
// 				// fmt.Println(val.Indel,val.IndelRef, val.IndelAlt, val.APos)
// 				// fmt.Println(idxVal, idxRes)
// 				switch val.Indel {
// 				case 0:
// 					PosInGenome[fname][val.APos] = fmt.Sprintf("%v/%v", val.NucInPosCoding, val.Alt)
// 					idxVal, idxRes := amino.GetComplexIndex(val.RefAAShort, val.AltAAShort, false)
// 					if val.TypeOf == "CDS" {
// 						// fmt.Printf("%v\t%v%v%v\t%v\t%v\t%v\t%v\t", val.Name, val.RefAA, val.CodonNbrInG, val.AltAA, val.Direction,
// 						// 	val.Mutation, idxVal, idxRes)
// 						// fmt.Println(idxVal)
// 						PosAdditionalInfo[val.APos] = fmt.Sprintf(
// 							"%v\t%v%v%v\t%v\t%v\t'%v\t%v\t%v\t", val.Name, val.RefAA, val.CodonNbrInG, val.AltAA,
// 							val.Direction,
// 							val.Mutation, idxVal, idxRes, val.PosInCodonG+1)
// 					} else {
// 						PosAdditionalInfo[val.APos] = fmt.Sprintf(
// 							"%v\t%v%v%v\t%v\t%v\t%v\t%v\t%v\t", val.Name, "-", "-", "-", val.Direction, "-", "-", "-", "-")
// 					}
// 				case 1:
// 					PosInGenome[fname][val.APos] = fmt.Sprintf("%v/%v", val.IndelRef, val.IndelAlt)
// 					PosAdditionalInfo[val.APos] = fmt.Sprintf(
// 						"%v\t%v%v%v\t%v\t%v\t'%v\t%v\t%v\t", val.Name, "-", "-", "-", val.Direction, "-", "-", "-", "-")
// 				}
//
// 				// AllPosUnsort = append(AllPosUnsort, val.APos)
//
// 			}
//
// 			// for _, allpos := range AllPos {
//
// 			// 	if pos[allpos] != "" {
// 			// 		posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v/%v", getNucFromGenomePos(allpos), pos[allpos]))
// 			// 		// fmt.Println(posFreq[allpos])
// 			// 		if len(group) != 0 {
// 			// 			fileGroup[allpos] = append(fileGroup[allpos], group[fname])
// 			// 			if len(label) != 0 {
// 			// 				fileLabel[allpos] = append(fileLabel[allpos], label[fname])
// 			// 			}
// 			// 			// groupArr = append(groupArr, group[fname])
// 			// 			// groupBody = fmt.Sprintf("\t%v\t", group[fname])
// 			// 		} else {
// 			// 			groupBody = "\t"
// 			// 		}
//
// 			// 	} else {
// 			// 		posFreq[allpos] = append(posFreq[allpos], ".")
//
// 			// 	}
//
// 			// }
//
// 		}
// 		// AllPos = unique(AllPosUnsort)
// 		// // allLocuses = removeStringDuplicates(allLocusUnsort)
// 		// sort.Ints(AllPos)
// 		// fmt.Println(PosInGenome)
//
// 		for _, gname := range genomes {
// 			for _, allpos := range AllPos {
// 				mut := PosInGenome[gname][allpos]
// 				if mut != "" {
// 					posFreq[allpos] = append(posFreq[allpos], mut)
// 					if len(group) != 0 {
// 						fileGroup[allpos] = append(fileGroup[allpos], group[strings.ToUpper(gname)])
//
// 						if len(label) != 0 {
// 							fileLabel[allpos] = append(fileLabel[allpos], label[strings.ToUpper(gname)])
// 						}
// 						// fmt.Printf("%v\t%v\t%v\n", allpos, group[strings.ToUpper(gname)], strings.ToUpper(gname))
//
// 					} else {
// 						groupBody = "\t"
// 					}
// 				} else {
// 					posFreq[allpos] = append(posFreq[allpos], ".")
// 				}
//
// 			}
//
// 		}
//
// 		for _, allpos := range AllPos {
// 			// if buffer.Len() == 0 {
// 			gname, _ := bsatseq.GetGeneNameByPos(allpos, allpos)
// 			if len(group) != 0 {
// 				sort.Strings(fileGroup[allpos])
// 				uniq := bsatservice.RmStrDoubles(fileGroup[allpos])
// 				if len(uniq) == 1 {
// 					if uniq[0] == "" {
// 						uniq = append(uniq, "-")
// 					}
// 				}
// 				groupBody = fmt.Sprintf("\t%v\t", strings.Trim(strings.Join(uniq, " "), " "))
//
// 			}
//
// 			if len(label) != 0 {
// 				sort.Strings(fileLabel[allpos])
// 				uniq := bsatservice.RmStrDoubles(fileLabel[allpos])
//
// 				groupBody = fmt.Sprintf("%v%v\t", groupBody, strings.Trim(strings.Join(uniq, " "), " "))
//
// 			}
//
// 			buffer.WriteString(fmt.Sprintln(allpos, "\t", gname, groupBody, PosAdditionalInfo[allpos], strings.Join(posFreq[allpos], "\t")))
// 			// fmt.Println(buffer.String())
//
// 		}
// 		headers.WriteString("\n")
// 		matrixPrint(headers, buffer, fileOut)
// 		// if buffer.Len() != 0 && headers.Len() != 0 {
// 		// 	fOut, err := os.Create(fileOut)
// 		// 	if err != nil {
// 		// 		log.Fatal("Cannot create file", err)
// 		// 	}
// 		// 	defer fOut.Close()
// 		// 	fmt.Fprintf(fOut, headers.String())
// 		// 	fmt.Fprintf(fOut, buffer.String())
// 		// 	fmt.Printf("\n\nWell done!\n")
// 		// t1 := time.Now()
// 		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
// 		// }
//
// 	case "locus":
// 		var (
// 			i          int
// 			locFreq    = map[string][]string{}
// 			locusCount = make(map[string]int)
// 		)
// 		headers.WriteString("Locus\t")
//
// 		// for i, file := range files {
// 		for fname, snps := range bsatstruct.SnpCacheMap {
// 			i++
// 			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))
//
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
//
// 			// snps := parserVCF(file, false, allGenesVal)
//
// 			for _, val := range snps {
//
// 				locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] = locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] + 1
//
// 			}
// 			for _, allloc := range allLocuses {
//
// 				if *statMinLocusCount != 0 {
// 					if locusCount[fmt.Sprintf("%v_%v", allloc, i)] >= *statMinLocusCount {
// 						locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))
// 					}
// 				} else {
// 					locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))
// 				}
//
// 			}
//
// 		}
//
// 		for _, allloc := range allLocuses {
//
// 			buffer.WriteString(fmt.Sprintln(allloc, "\t", strings.Join(locFreq[allloc], "\t")))
//
// 		}
//
// 		headers.WriteString("\n")
// 		matrixPrint(headers, buffer, fileOut)
//
// 	case "freq":
// 		headers.WriteString("Pos\tFreq\n")
// 		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
// 		for _, snps := range bsatstruct.SnpCacheMap {
// 			for _, val := range snps {
//
// 				posCount[val.APos] = posCount[val.APos] + 1
//
// 			}
//
// 		}
//
// 		for _, allpos := range AllPos {
//
// 			// headers.WriteString(fmt.Sprintf("P%v\n", allpos))
// 			buffer.WriteString(fmt.Sprintf("%v\t%v\n", allpos, posCount[allpos]))
// 			// fmt.Println(posCount[allpos])
// 		}
//
// 		matrixPrint(headers, buffer, fileOut)
//
// 	case "dnds":
// 		matrixDnDs(fileOut)
// 	case "mst":
// 		// stat --db test_core -a matrix -t mst -o test.tsv --exclude-genes=exgenes.txt --exclude-snp=drugs2.txt --snp-number=15
// 		var (
// 			exGenes = make(map[int]int)
// 			exSNP   = make(map[int]int)
// 		)
// 		if *gbExcludeGenes != "" {
// 			exGenes = bsatf.LoadExclGenes(*gbExcludeGenes)
// 		} else if *gbExcludeSnp != "" {
// 			exSNP = bsatf.LoadExclSNP(*gbExcludeSnp)
// 		}
//
// 		matrixBinaryMST(fileOut, *gbVerbose, exGenes, exSNP, *gbRandomize)
//
// 	case "jw":
// 		// var dnds [][]DnDsRes
// 		var (
// 			locJW = map[string][]string{}
// 			i     int
// 		)
// 		headers.WriteString("Locus\t")
//
// 		// for i, file := range files {
// 		for fname := range bsatstruct.SnpCacheMap {
// 			i++
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
//
// 			fmt.Printf("Calculating Jaro Winkler distance: Working on %v from %v (%v) \r", i, len(*files), fname)
//
// 			jw := bsatcalc.CalcJaroWinklerDist(fname, *gbVerbose)
//
// 			for _, jwVal := range jw {
//
// 				locJW[jwVal.Locus] = append(locJW[jwVal.Locus], fmt.Sprintf("%.3f", jwVal.JWDist))
//
// 			}
//
// 		}
//
// 		for _, allloc := range allLocuses {
//
// 			// if len(locJW[allloc]) != 0 {
//
// 			// fmt.Println(len(*files)-len(locJW[allloc]), locJW[allloc], strings.Repeat(" 1 ", len(*files)-len(locJW[allloc])))
// 			rpt := strings.Repeat("1\t", len(*files)-len(locJW[allloc]))
// 			locJW[allloc] = append(locJW[allloc], rpt)
// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t", allloc), strings.Join(locJW[allloc], "\t")))
//
// 			// }
// 		}
//
// 		headers.WriteString("\n")
// 		matrixPrint(headers, buffer, fileOut)
//
// 	case "gc3":
// 		// calcGC3Val
//
// 		var (
// 			locGC3 = map[string][]string{}
// 			refGC3 = map[string]string{}
// 			i      int
// 		)
// 		headers.WriteString("Locus\tRefCG3\t")
//
// 		// for i, file := range files {
// 		for fname := range bsatstruct.SnpCacheMap {
// 			// fmt.Println(fname)
// 			i++
// 			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
//
// 			fmt.Printf("Calculating GC3 values: Working on %v from %v\t%v \r", i, len(*files), fname)
// 			// logger.Printf("%v\n", fname)
// 			// fmt.Println(snpCacheMap[fname])
// 			gc3 := bsatcalc.CalcGC3Val(bsatstruct.SnpCacheMap[fname])
// 			// logger.Printf("%v\t%v\n", fname, gc3)
// 			// fmt.Println(gc3)
// 			// fmt.Println(gc3)
//
// 			for _, gc3val := range gc3 {
//
// 				locGC3[gc3val.Locus] = append(locGC3[gc3val.Locus], fmt.Sprintf("%.2f", gc3val.GC3Alt))
// 				refGC3[gc3val.Locus] = fmt.Sprintf("%.2f", gc3val.GC3Ref)
//
// 				// if *flgDebug == true {
// 				// 	fmt.Printf("L:%v Alt:%v Ref:%v\n", gc3val.Locus, gc3val.GC3Alt, gc3val.GC3Ref)
// 				// }
//
// 			}
//
// 		}
//
// 		for _, allloc := range allLocuses {
// 			// if len(locGC3[allloc]) != 0 {
// 			if len(locGC3[allloc]) != 0 {
// 				rpt := strings.Repeat(refGC3[allloc]+"\t", len(*files)-len(locGC3[allloc]))
// 				locGC3[allloc] = append(locGC3[allloc], rpt)
// 				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC3[allloc]), strings.Join(locGC3[allloc], "\t")))
//
// 			}
// 		}
//
// 		headers.WriteString("\n")
// 		matrixPrint(headers, buffer, fileOut)
//
// 	}
//
// }

func MakeMatrix(typeof string, fileOut string, verbose bool) {

	var (
		AllPos     []int
		allLocuses []string
		buffer     strings.Builder
		headers    strings.Builder
		posCount   = make(map[int]int)
	)

	files := &bsatstruct.ListOfFiles

	for i, file := range *files {
		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery)}
		go qSNP.Request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for key, val := range bsatstruct.GenePositions {
		if val.Type == "CDS" {
			allLocuses = append(allLocuses, key)
		}
	}
	sort.Strings(allLocuses)

	var (
		exGenes = make(map[int]int)
		exSNP   = make(map[int]int)
	)
	if bsatstruct.Flag.GbExcludeGenes != "" {
		exGenes = bsatf.LoadExclGenes(bsatstruct.Flag.GbExcludeGenes)
	} else if bsatstruct.Flag.GbExcludeSnp != "" {
		exSNP = bsatf.LoadExclSNP(bsatstruct.Flag.GbExcludeSnp)
	}

	switch typeof {

	case "binary":

		MatrixBinary(fileOut, exGenes, exSNP)

	// case "table":

	// 	matrixTable(fileOut)

	case "nc":
		var posFreq = map[int][]string{}
		pos := make(map[int]string)

		headers.WriteString("Pos\tRef\t")
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
		for fname, snps := range bsatstruct.SnpCacheMap {
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt
				// AllPosUnsort = append(AllPosUnsort, val.APos)

			}
			for _, allpos := range AllPos {

				if pos[allpos] != "" {
					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", pos[allpos]))
					// fmt.Println(posFreq[allpos])

				} else {
					posFreq[allpos] = append(posFreq[allpos], bsatseq.GetNucFromGenomePos(allpos))

				}

			}

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			buffer.WriteString(fmt.Sprintln(allpos, "\t", bsatseq.GetNucFromGenomePos(allpos), "\t", strings.Join(posFreq[allpos], "\t")))
			// fmt.Println(buffer.String())

		}
		headers.WriteString("\n")
		bsatservice.MatrixPrint(headers, buffer, fileOut)
	case "nc_coded":
		// A, T, G,C, ---> 0, 1, 2, 3
		var posFreq = map[int][]string{}
		pos := make(map[int]string)
		posCount := make(map[int]int)
		files := &bsatstruct.ListOfFiles
		headers.WriteString("Pos\t")
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
		for fname, snps := range bsatstruct.SnpCacheMap {
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				pos[val.APos] = val.Alt
				// AllPosUnsort = append(AllPosUnsort, val.APos)
				posCount[val.APos] = posCount[val.APos] + 1
			}
			for _, allpos := range AllPos {

				if pos[allpos] != "" {

					posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v", bsatseq.Nuc2IntCode(pos[allpos])))
					// fmt.Println(posFreq[allpos])

				} else {
					posFreq[allpos] = append(posFreq[allpos], bsatseq.Nuc2IntCode(bsatseq.GetNucFromGenomePos(allpos)))

				}

			}

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			if posCount[allpos] <= len(*files)-1 {
				buffer.WriteString(fmt.Sprintln(allpos, "\t", strings.Join(posFreq[allpos], "\t")))
				// fmt.Println(buffer.String())
			}

		}
		headers.WriteString("\n")
		bsatservice.MatrixPrint(headers, buffer, fileOut)

	case "summary":
		var (
			posFreq = map[int][]string{}
			// pos                    = make(map[int]string)
			groupRegexp            = regexp.MustCompile(`^(\S+)\W+(\w+)\W+(\w+)`)
			group                  = make(map[string]string)
			label                  = make(map[string]string)
			groupHeader, groupBody string
			fileGroup              = make(map[int][]string)
			fileLabel              = make(map[int][]string)
			PosInGenome            = make(map[string]map[int]string)
			PosAdditionalInfo      = make(map[int]string)
			genomes                []string
		)

		if bsatstruct.Flag.StatGroupFromFile != "" {
			f, err := os.Open(bsatstruct.Flag.StatGroupFromFile) // открываем файл

			if err != nil {
				log.Fatal(err)

			}
			defer func(f *os.File) {
				err := f.Close()
				if err != nil {
					log.Fatal(err)
				}
			}(f)

			scanner := bufio.NewScanner(f) //  новый сканер

			for scanner.Scan() {

				scanTxt := scanner.Text()

				for _, grpVal := range groupRegexp.FindAllStringSubmatch(scanTxt, -1) {
					group[strings.ToUpper(grpVal[1])] = strings.ToUpper(grpVal[2])

					label[strings.ToUpper(grpVal[1])] = grpVal[3]

					// fmt.Println(grpVal)
				}
			}

		}

		// fmt.Println(group)
		if len(group) != 0 {
			groupHeader = "Group\t"
			if len(label) != 0 {
				groupHeader = fmt.Sprintf("%vLabel\t", groupHeader)
			}
		}

		headers.WriteString(fmt.Sprintf("Pos\tLocus\t%vGene\tAAchange\tDirection\tMutation\tCIndex\tCIndex_Res\tPosInCodon\t", groupHeader))
		// allPosChan := make(chan []int)
		// go func() {
		// 	allPosChan <- getAllPosFromCacheMap()
		// }()
		// AllPos = <-allPosChan
		// for i, file := range files {
		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)

		for fname, snps := range bsatstruct.SnpCacheMap {
			PosInGenome[fname] = make(map[int]string)
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))
			genomes = append(genomes, fname)
			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				// pos[val.APos] = val.Alt

				// fmt.Println(val.Indel,val.IndelRef, val.IndelAlt, val.APos)
				// fmt.Println(idxVal, idxRes)
				switch val.Indel {
				case 0:
					PosInGenome[fname][val.APos] = fmt.Sprintf("%v/%v", val.NucInPosCoding, val.Alt)
					idxVal, idxRes := amino.GetComplexIndex(val.RefAAShort, val.AltAAShort, false)
					if val.TypeOf == "CDS" {
						// fmt.Printf("%v\t%v%v%v\t%v\t%v\t%v\t%v\t", val.Name, val.RefAA, val.CodonNbrInG, val.AltAA, val.Direction,
						// 	val.Mutation, idxVal, idxRes)
						// fmt.Println(idxVal)
						PosAdditionalInfo[val.APos] = fmt.Sprintf(
							"%v\t%v%v%v\t%v\t%v\t'%v\t%v\t%v\t", val.Name, val.RefAA, val.CodonNbrInG, val.AltAA,
							val.Direction,
							val.Mutation, idxVal, idxRes, val.PosInCodonG+1)
					} else {
						PosAdditionalInfo[val.APos] = fmt.Sprintf(
							"%v\t%v%v%v\t%v\t%v\t%v\t%v\t%v\t", val.Name, "-", "-", "-", val.Direction, "-", "-", "-", "-")
					}
				case 1:
					PosInGenome[fname][val.APos] = fmt.Sprintf("%v/%v", val.IndelRef, val.IndelAlt)
					PosAdditionalInfo[val.APos] = fmt.Sprintf(
						"%v\t%v%v%v\t%v\t%v\t'%v\t%v\t%v\t", val.Name, "-", "-", "-", val.Direction, "-", "-", "-", "-")
				}

				// AllPosUnsort = append(AllPosUnsort, val.APos)

			}

			// for _, allpos := range AllPos {

			// 	if pos[allpos] != "" {
			// 		posFreq[allpos] = append(posFreq[allpos], fmt.Sprintf("%v/%v", getNucFromGenomePos(allpos), pos[allpos]))
			// 		// fmt.Println(posFreq[allpos])
			// 		if len(group) != 0 {
			// 			fileGroup[allpos] = append(fileGroup[allpos], group[fname])
			// 			if len(label) != 0 {
			// 				fileLabel[allpos] = append(fileLabel[allpos], label[fname])
			// 			}
			// 			// groupArr = append(groupArr, group[fname])
			// 			// groupBody = fmt.Sprintf("\t%v\t", group[fname])
			// 		} else {
			// 			groupBody = "\t"
			// 		}

			// 	} else {
			// 		posFreq[allpos] = append(posFreq[allpos], ".")

			// 	}

			// }

		}
		// AllPos = unique(AllPosUnsort)
		// // allLocuses = removeStringDuplicates(allLocusUnsort)
		// sort.Ints(AllPos)
		// fmt.Println(PosInGenome)

		for _, gname := range genomes {
			for _, allpos := range AllPos {
				mut := PosInGenome[gname][allpos]
				if mut != "" {
					posFreq[allpos] = append(posFreq[allpos], mut)
					if len(group) != 0 {
						fileGroup[allpos] = append(fileGroup[allpos], group[strings.ToUpper(gname)])

						if len(label) != 0 {
							fileLabel[allpos] = append(fileLabel[allpos], label[strings.ToUpper(gname)])
						}
						// fmt.Printf("%v\t%v\t%v\n", allpos, group[strings.ToUpper(gname)], strings.ToUpper(gname))

					} else {
						groupBody = "\t"
					}
				} else {
					posFreq[allpos] = append(posFreq[allpos], ".")
				}

			}

		}

		for _, allpos := range AllPos {
			// if buffer.Len() == 0 {
			gname, _ := bsatseq.GetGeneNameByPos(allpos, allpos)
			if len(group) != 0 {
				sort.Strings(fileGroup[allpos])
				uniq := bsatservice.RmStrDoubles(fileGroup[allpos])
				if len(uniq) == 1 {
					if uniq[0] == "" {
						uniq = append(uniq, "-")
					}
				}
				groupBody = fmt.Sprintf("\t%v\t", strings.Trim(strings.Join(uniq, " "), " "))

			}

			if len(label) != 0 {
				sort.Strings(fileLabel[allpos])
				uniq := bsatservice.RmStrDoubles(fileLabel[allpos])

				groupBody = fmt.Sprintf("%v%v\t", groupBody, strings.Trim(strings.Join(uniq, " "), " "))

			}

			buffer.WriteString(fmt.Sprintln(allpos, "\t", gname, groupBody, PosAdditionalInfo[allpos], strings.Join(posFreq[allpos], "\t")))
			// fmt.Println(buffer.String())

		}
		headers.WriteString("\n")
		bsatservice.MatrixPrint(headers, buffer, fileOut)
		// if buffer.Len() != 0 && headers.Len() != 0 {
		// 	fOut, err := os.Create(fileOut)
		// 	if err != nil {
		// 		log.Fatal("Cannot create file", err)
		// 	}
		// 	defer fOut.Close()
		// 	fmt.Fprintf(fOut, headers.String())
		// 	fmt.Fprintf(fOut, buffer.String())
		// 	fmt.Printf("\n\nWell done!\n")
		// t1 := time.Now()
		// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
		// }

	case "locus":
		var (
			i          int
			locFreq    = map[string][]string{}
			locusCount = make(map[string]int)
		)
		headers.WriteString("Locus\t")

		// for i, file := range files {
		for fname, snps := range bsatstruct.SnpCacheMap {
			i++
			// fmt.Printf("Generating matrix: Working on  %v from %v \r", i+1, len(files))

			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			// snps := parserVCF(file, false, allGenesVal)

			for _, val := range snps {

				locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] = locusCount[fmt.Sprintf("%v_%v", val.Locus, i)] + 1

			}
			for _, allloc := range allLocuses {

				if bsatstruct.Flag.StatMinLocusCount != 0 {
					if locusCount[fmt.Sprintf("%v_%v", allloc, i)] >= bsatstruct.Flag.StatMinLocusCount {
						locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))
					}
				} else {
					locFreq[allloc] = append(locFreq[allloc], strconv.Itoa(locusCount[fmt.Sprintf("%v_%v", allloc, i)]))
				}

			}

		}

		for _, allloc := range allLocuses {

			buffer.WriteString(fmt.Sprintln(allloc, "\t", strings.Join(locFreq[allloc], "\t")))

		}

		headers.WriteString("\n")
		bsatservice.MatrixPrint(headers, buffer, fileOut)

	case "freq":
		headers.WriteString("Pos\tFreq\n")
		AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
		for _, snps := range bsatstruct.SnpCacheMap {
			for _, val := range snps {

				posCount[val.APos] = posCount[val.APos] + 1

			}

		}

		for _, allpos := range AllPos {

			// headers.WriteString(fmt.Sprintf("P%v\n", allpos))
			buffer.WriteString(fmt.Sprintf("%v\t%v\n", allpos, posCount[allpos]))
			// fmt.Println(posCount[allpos])
		}

		bsatservice.MatrixPrint(headers, buffer, fileOut)

	case "dnds":
		MatrixDnDs(fileOut)
	case "mst":
		// stat --db test_core -a matrix -t mst -o test.tsv --exclude-genes=exgenes.txt --exclude-snp=drugs2.txt --snp-number=15
		var (
			exGenes = make(map[int]int)
			exSNP   = make(map[int]int)
		)
		if bsatstruct.Flag.GbExcludeGenes != "" {
			exGenes = bsatf.LoadExclGenes(bsatstruct.Flag.GbExcludeGenes)
		} else if bsatstruct.Flag.GbExcludeSnp != "" {
			exSNP = bsatf.LoadExclSNP(bsatstruct.Flag.GbExcludeSnp)
		}

		MatrixBinaryMST(fileOut, bsatstruct.Flag.GbVerbose, exGenes, exSNP, bsatstruct.Flag.GbRandomize)

	case "jw":
		// var dnds [][]DnDsRes
		var (
			locJW = map[string][]string{}
			i     int
		)
		headers.WriteString("Locus\t")

		// for i, file := range files {
		for fname := range bsatstruct.SnpCacheMap {
			i++
			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			fmt.Printf("Calculating Jaro Winkler distance: Working on %v from %v (%v) \r", i, len(*files), fname)

			jw := CalcJaroWinklerDist(fname, bsatstruct.Flag.GbVerbose)

			for _, jwVal := range jw {

				locJW[jwVal.Locus] = append(locJW[jwVal.Locus], fmt.Sprintf("%.3f", jwVal.JWDist))

			}

		}

		for _, allloc := range allLocuses {

			// if len(locJW[allloc]) != 0 {

			// fmt.Println(len(*files)-len(locJW[allloc]), locJW[allloc], strings.Repeat(" 1 ", len(*files)-len(locJW[allloc])))
			rpt := strings.Repeat("1\t", len(*files)-len(locJW[allloc]))
			locJW[allloc] = append(locJW[allloc], rpt)
			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t", allloc), strings.Join(locJW[allloc], "\t")))

			// }
		}

		headers.WriteString("\n")
		bsatservice.MatrixPrint(headers, buffer, fileOut)

	case "gc3":
		// calcGC3Val

		var (
			locGC3 = map[string][]string{}
			refGC3 = map[string]string{}
			i      int
		)
		headers.WriteString("Locus\tRefCG3\t")

		// for i, file := range files {
		for fname := range bsatstruct.SnpCacheMap {
			// fmt.Println(fname)
			i++
			headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(fname, filepath.Ext(fname))))

			fmt.Printf("Calculating GC3 values: Working on %v from %v\t%v \r", i, len(*files), fname)
			// logger.Printf("%v\n", fname)
			// fmt.Println(snpCacheMap[fname])
			gc3 := CalcGC3Val(bsatstruct.SnpCacheMap[fname])
			// logger.Printf("%v\t%v\n", fname, gc3)
			// fmt.Println(gc3)
			// fmt.Println(gc3)

			for _, gc3val := range gc3 {

				locGC3[gc3val.Locus] = append(locGC3[gc3val.Locus], fmt.Sprintf("%.2f", gc3val.GC3Alt))
				refGC3[gc3val.Locus] = fmt.Sprintf("%.2f", gc3val.GC3Ref)

				// if *flgDebug == true {
				// 	fmt.Printf("L:%v Alt:%v Ref:%v\n", gc3val.Locus, gc3val.GC3Alt, gc3val.GC3Ref)
				// }

			}

		}

		for _, allloc := range allLocuses {
			// if len(locGC3[allloc]) != 0 {
			if len(locGC3[allloc]) != 0 {
				rpt := strings.Repeat(refGC3[allloc]+"\t", len(*files)-len(locGC3[allloc]))
				locGC3[allloc] = append(locGC3[allloc], rpt)
				buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v\t%v\t", allloc, refGC3[allloc]), strings.Join(locGC3[allloc], "\t")))

			}
		}

		headers.WriteString("\n")
		bsatservice.MatrixPrint(headers, buffer, fileOut)

	}

}

func MatrixBinary(fileOut string, exGenes map[int]int, exSNP map[int]int) {

	/*
		создание бинарной матрицы

	*/

	var (
		AllPosUnsort, AllPos []int
		allLocusUnsort       []string
		buffer               strings.Builder
		headers              strings.Builder
		posCount             = make(map[int]int)
		// snps                 []snpInfo
		posFN   = make(map[int][]string)
		posFreq = map[int][]string{}
		// resPOS  = make(map[int]int)
		th int
		// countPassSNPfromSNPList, countPassSNPinGenes int
	)

	// type binMatrixStruct struct {
	// 	pos, res int
	// 	posArray []string
	// }
	// var ResSeq []seqInfo

	// files := GetListVcf()
	files := &bsatstruct.ListOfFiles

	if bsatstruct.Flag.StatTH < len(*files) {
		th = bsatstruct.Flag.StatTH
		// fmt.Println("!", th)
	} else if bsatstruct.Flag.StatTH >= len(*files) {
		th = len(*files) - 1
		// fmt.Println("!!", th)

	}

	// fmt.Println(files)
	pos := make(map[int]string)

	headers.WriteString("\t")
	// i := 1
	for key, snps := range bsatstruct.SnpCacheMap {
		headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))
		//
		// fmt.Printf("Counting SNP positions: Working on %v files from %v \r", i+1, len(*files))

		for _, val := range snps {

			pos[val.APos] = val.Alt
			// posTest = append(posTest, posByFile{pos: pos, file: file, apos: val.APos})
			if strings.Contains(strings.Join(posFN[val.APos], " "), key) == false {
				posFN[val.APos] = append(posFN[val.APos], key)

			}

			AllPosUnsort = append(AllPosUnsort, val.APos)
			posCount[val.APos] = posCount[val.APos] + 1
			if val.TypeOf == "CDS" {
				allLocusUnsort = append(allLocusUnsort, val.Locus)
			}
		}

	}
	// fmt.Println("Passed SNP from Genes: ", countPassSNPinGenes, " , passed from SNP list: ", countPassSNPfromSNPList)
	AllPos = bsatseq.GetAllPositions(exGenes, exSNP)
	// allLocuses = removeStringDuplicates(allLocusUnsort)
	sort.Ints(AllPos)

	for _, file := range *files {
		for _, allpos := range AllPos {

			if strings.Contains(strings.Join(posFN[allpos], " "), file) {
				posFreq[allpos] = append(posFreq[allpos], "1")

			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}

	for _, pos := range AllPos {
		var count0, count1 int
		for i := 0; i < len(posFreq[pos]); i++ {
			if posFreq[pos][i] == "0" {
				count0++
			} else if posFreq[pos][i] == "1" {
				count1++
			}
		}
		// resPOS[pos] = n
		if bsatstruct.Flag.StatTH == 0 {
			// if *gbDebug {
			// 	fmt.Println(count0, count1, th, posFreq[pos], len(posFreq[pos]))
			// }
			if count1 <= len(posFreq[pos])-1 && count0 <= len(posFreq[pos])-1 {
				buffer.WriteString(fmt.Sprintln(pos, "\t", strings.Join(posFreq[pos], "\t")))
			}
		} else {

			if count0 >= th && count1 >= th && th <= len(posFreq[pos]) && count1 < len(posFreq[pos]) && count1 < len(posFreq[pos]) {
				if bsatstruct.Flag.GbDebug == true {
					fmt.Printf("pos:%v isRef:%v,isAlt:%v,th:%v %v legnth_array: %v\n", pos, count0, count1, th, posFreq[pos], len(posFreq[pos]))
				}
				buffer.WriteString(fmt.Sprintln(pos, "\t", strings.Join(posFreq[pos], "\t")))
			}

		}

	}

	headers.WriteString("\n")
	bsatservice.MatrixPrint(headers, buffer, fileOut)

}

func MatrixBinaryMST(fileOut string, verbose bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) {

	var (
		AllPos, SelectedPos []int
		// ResSeq              []seqInfo
		// passSNP = make(map[string]int)
		uniqSNP         = make(map[int]int)
		nbrOfSNP        = 5000
		headers, buffer strings.Builder
		snpCount        = make(map[int]int)
		// posCount             = make(map[int]int)
	)
	// files := GetListVcf()

	// queryChan := make(chan vcfInfoQuery)

	if bsatstruct.Flag.StatNbrOfSNP != 0 {
		nbrOfSNP = bsatstruct.Flag.StatNbrOfSNP
	}

	files := &bsatstruct.ListOfFiles
	for i, file := range *files {
		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery), Print: verbose}
		go qSNP.Request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for _, snps := range bsatstruct.SnpCacheMap {

		for _, val := range snps {
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// fmt.Println(val.Locus, val.Start, val.End, val.Product)
						// passSNP["genes"] = passSNP["genes"] + 1
						uniqSNP[val.APos] = 2 // 2-EXCLUDED
						continue
					} else if exSNPs[val.APos] == 1 {
						// passSNP["snp"] = passSNP["snp"] + 1
						uniqSNP[val.APos] = 2 // 2-EXCLUDED
						continue

						// fmt.Println(val.APos)

					} else {
						if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
							uniqSNP[val.APos] = 1 // 1-INCLUDED

						}
					}
				}
			} else {
				if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
					uniqSNP[val.APos] = 1

				}
				// AllPosUnsort = append(AllPosUnsort, val.APos)
			}
			snpCount[val.APos] = snpCount[val.APos] + 1
			// if exGenes[val.Locus] != 1 && exSNPs[string(val.APos)] != 1 {
			// 	AllPosUnsort = append(AllPosUnsort, val.APos)
			// }
		}
	}

	// fmt.Println(uniqSNP)

	for key, value := range uniqSNP {

		if value == 1 && len(*files) != snpCount[key] {

			AllPos = append(AllPos, key)
			// } else if value == 2 {
			// 	fmt.Println(key)

		}

	}
	// go process("Working...          ")
	// AllPos = unique(AllPosUnsort)

	sort.Ints(AllPos)

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true {
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[rnd])
			}

		}
	} else {
		for i := 0; i <= nbrOfSNP; i++ {
			// rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[i])
			}

		}
	}
	// fmt.Println(AllPos[0])

	// rand.Intn(max - min) + min

	sort.Ints(SelectedPos)

	// for i := 0; i < len(SelectedPos); i++ {
	// 	headers.WriteString(fmt.Sprintln(SelectedPos[i], "\t"))
	// }
	// headers.WriteString(fmt.Sprintf("%v\n", strings.Join(SelectedPos, "\t")))
	headers.WriteString("FILE\t")
	// buffer.WriteString("\t")
	headers.WriteString(strings.Trim(strings.Join(strings.Fields(fmt.Sprint(SelectedPos)), "\t"), "[]"))
	headers.WriteString("\n")
	if verbose == true {
		fmt.Println(headers.String())
		// fmt.Println(AllPos)
	}
	for fname, snps := range bsatstruct.SnpCacheMap {
		buffer.WriteString(fname)
		// if verbose == true {
		// 	fmt.Printf("Generating sequences: Working on  %v from %v \r", i+1, len(snps.File))
		// }
		pos := make(map[int]string)
		// var buffer strings.Builder

		// buffer.WriteString(fmt.Sprintf(">%v\n", strings.ToUpper(fname)))

		for _, val := range snps {
			pos[val.APos] = val.Alt
			// if exGenes[val.Locus] == 1 {
			// 	fmt.Println(val.Locus)
			// }

		}
		for _, allpos := range SelectedPos {
			// posCount[allpos] = posCount[allpos] + 1
			// fmt.Println(snpCount[allpos])

			if pos[allpos] != "" {

				buffer.WriteString("\t1")
			} else {

				buffer.WriteString("\t0")
			}
			// }

		}
		buffer.WriteString("\n")
		// ResSeq = append(ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: SelectedPos})
		// fmt.Println(len(buffer.String()))

		// if *gbDebug == true {
		// 	fmt.Printf("%v\t:\nThere was passed %v SNPs from exclude gene file\n And %v SNPs from exclude snp file\n", fname, passSNP["genes"], passSNP["snp"])
		// }
	}
	if verbose == true {
		fmt.Println(buffer.String())
	}
	bsatservice.MatrixPrint(headers, buffer, fileOut)
	// return ResSeq
}

func MatrixDnDs(fileOut string) {
	// var AllPos []int

	// type dndsPerGenome struct {
	// 	dnds, locus string
	// }
	var (
		allLocuses []string

		buffer  strings.Builder
		headers strings.Builder
		// var posCount = make(map[int]int)
		// snps         []snpInfo
		// altPositions        = make(map[string][]allPositionsInGene)
		altPositionsPerFile = make(map[string]map[string][]bsatstruct.AllPositionsInGene)
		// locDNDS                        = map[string][]string{}
		countNbrOne        = make(map[string]int)
		countNonZeroValues = make(map[string]int)
		countPositive      = make(map[string]int)
		countNegative      = make(map[string]int)
		countNeutral       = make(map[string]int)
		positiveGenes      = make(map[string][]string)
		geneDnDs           = make(map[string]map[string]string)
		positiveGenesList  []string
		positiveGenesCheck = make(map[string]int)

		// locInGenome                    = make(map[string]map[string]string)
		// genomes                        []string
		usedLocusesUnsort, usedLocuses []string
		dndsPerLocus                   = make(map[string][]string)
		// filePerLocus                   = make(map[string][]string)
		groupRegexp            = regexp.MustCompile(`^(\S+)\W+(\w+)\W+(\w+)`)
		group                  = make(map[string]string)
		label                  = make(map[string]string)
		groupHeader, groupBody string
		fileGroup              = make(map[string][]string)
		fileLabel              = make(map[string][]string)
		dndsPerGenomeAll       = make(map[string]map[string]string)
		t0                     = time.Now()
		th                     int
		dndsFloat              float64
	)

	if bsatstruct.Flag.StatDnDsCountTh != 0 {
		th = bsatstruct.Flag.StatDnDsCountTh
	} else {
		th = 1
	}
	files := &bsatstruct.ListOfFiles
	// fmt.Println(files)

	i := 1
	// headers.WriteString("Locus\t")

	if bsatstruct.Flag.StatGroupFromFile != "" {
		f, err := os.Open(bsatstruct.Flag.StatGroupFromFile) // открываем файл

		if err != nil {
			log.Fatal(err)

		}
		defer func(f *os.File) {
			err := f.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(f)

		scanner := bufio.NewScanner(f) //  новый сканер

		for scanner.Scan() {

			scanTxt := scanner.Text()

			for _, grpVal := range groupRegexp.FindAllStringSubmatch(scanTxt, -1) {
				group[grpVal[1]] = grpVal[2]
				label[grpVal[1]] = grpVal[3]

				// fmt.Println(grpVal)
			}
		}

	}
	// fmt.Println(group)
	if len(group) != 0 {
		groupHeader = "Group\t"
		if len(label) != 0 {
			groupHeader = fmt.Sprintf("%vLabel\tProduct\t", groupHeader)
		}

	} else {

		groupHeader = "Product\t"
	}

	headers.WriteString(fmt.Sprintf("Locus\t%v", groupHeader))

	for key, val := range bsatstruct.GenePositions {
		if val.Type == "CDS" {
			allLocuses = append(allLocuses, key)
		}

	}
	// fmt.Println(geneCoordinates)
	for fname, snps := range bsatstruct.SnpCacheMap {
		// t0 = time.Now()
		// locInGenome[fname] = make(map[string]string)
		// headers.WriteString(fmt.Sprintf("%v\t", strings.TrimSuffix(key, filepath.Ext(key))))

		// altPositions = nil

		// genomes = append(genomes, fname)

		fmt.Printf("Gathering information: Working on %v from %v (%v)\r", i, len(*files), fname)
		altPositionsPerFile[fname] = make(map[string][]bsatstruct.AllPositionsInGene)
		for _, val := range snps {

			if val.TypeOf == "CDS" && bsatseq.ChkPosExists(altPositionsPerFile[fname][val.Locus], val.PosInGene, val.Alt) == false {

				// altPositions[val.Locus] = append(altPositions[val.Locus], allPositionsInGene{pos: val.PosInGene, alt: val.Alt, ref: val.NucInPosCoding, locus: val.Locus})
				altPositionsPerFile[fname][val.Locus] = append(
					altPositionsPerFile[fname][val.Locus], bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding,
						Locus: val.Locus})
				usedLocusesUnsort = append(usedLocusesUnsort, val.Locus)
			}

		}
	}
	// fmt.Println(altPositionsPerFile)
	// dndsPerGenomeAll = make(map[string]string)
	usedLocuses = bsatservice.RmStrDoubles(usedLocusesUnsort)
	for _, fname := range *files {
		dndsPerGenomeAll[fname] = make(map[string]string)
		geneDnDs[fname] = map[string]string{}
		t1 := time.Now()
		for _, allloc := range allLocuses {

			// prod := getProductByName(allloc)
			if len(altPositionsPerFile[fname][allloc]) > 2 {
				dndsChan := make(chan []string)
				go func() {

					dndsChan <- GetDnDsByLocus(allloc, altPositionsPerFile[fname][allloc])
				}()
				dndsRes, ok := <-dndsChan

				// fmt.Println(fname, dndsRes)
				if ok {

					// if allloc == "Rv3879c" {

					// 	fmt.Println(allloc, dndsRes, altPositions[allloc], fname)

					// }

					// fmt.Println(dndsRes)

					// locDNDS[dndsRes[0]] = append(locDNDS[dndsRes[0]], dndsRes[1])
					// if allloc == "Rv3854c" {
					// 	fmt.Println(allloc, dndsRes, altPositions[allloc])
					// }
					// fmt.Println(allloc, dndsRes)
					// locInGenome[fname][dndsRes[0]] = dndsRes[1]
					if dndsRes[1] == "1.00" {
						countNbrOne[allloc]++
						countNeutral[fname]++
					}
					// close(dndsChan)
					dndsPerGenomeAll[fname][dndsRes[0]] = dndsRes[1]
					dndsFloat, _ = strconv.ParseFloat(dndsRes[1], 64)
					geneDnDs[fname][allloc] = dndsRes[1]
					if dndsFloat > 1 {
						countPositive[fname]++
						positiveGenes[fname] = append(positiveGenes[fname], allloc)
						positiveGenesList = append(positiveGenesList, allloc)
						positiveGenesCheck[allloc] = 1
					} else if dndsFloat < 1 {
						countNegative[fname]++
					}
					countNonZeroValues[allloc]++
				}

			} else {

				// locDNDS[allloc] = append(locDNDS[allloc], "1")
				// locInGenome[fname][allloc] = "1.00"
				dndsPerGenomeAll[fname][allloc] = "1.00"
				geneDnDs[fname][allloc] = "1.00"
				// geneDnDs[fname][allloc] = "1.00"
				countNeutral[fname]++
				countNbrOne[allloc]++
			}

		}

		// if *gbVerbose == true {
		fmt.Printf("Calculating DN/DS: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), fname, t1.Sub(t0))
		// }
		i++
		// fmt.Println(locInGenome)
		// fmt.Println(len(altPositions))

	}
	i = 1

	// sort.Strings(genomes)
	if bsatstruct.Flag.GbDebug == true {
		fmt.Println("Used locuses:", usedLocuses, "Total:", len(usedLocuses))
	}

	for _, gname := range *files {
		t1 := time.Now()
		for i := 0; i < len(usedLocuses); i++ {
			// for j := 0; j < len(dndsPerGenomeAll); j++ {
			// if dndsPerGenomeAll[gname] == usedLocuses[i] {
			dndsPerLocus[usedLocuses[i]] = append(dndsPerLocus[usedLocuses[i]], dndsPerGenomeAll[gname][usedLocuses[i]])
			// filePerLocus[allloc] = append(filePerLocus[allloc], gname)

			// }
			// }
		}
		// for _, allloc := range usedLocuses {
		// 	for _, dnds := range dndsPerGenomeAll {
		// 		if dnds.genome == gname && dnds.locus == allloc {
		// 			dndsPerLocus[allloc] = append(dndsPerLocus[allloc], dnds.dnds)
		// 			// filePerLocus[allloc] = append(filePerLocus[allloc], gname)

		// 		}
		// 	}
		// }

		fmt.Printf("Making matrix: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), gname, t1.Sub(t0))
		i++
	}

	i = 1

	for _, gname := range *files {
		// t1 := time.Now()
		headers.WriteString(fmt.Sprintf("%v\t", gname))
		for _, allloc := range usedLocuses {
			// dndsPerLocus[allloc] = append(dndsPerLocus[allloc], locInGenome[gname][allloc])

			if len(group) != 0 {
				sort.Strings(fileGroup[allloc])
				uniq := bsatservice.RmStrDoubles(fileGroup[allloc])
				groupBody = fmt.Sprintf("\t%v\t", strings.Join(uniq, ""))
				// fmt.Println(groupBody)
			}
			if len(label) != 0 {
				sort.Strings(fileLabel[allloc])
				uniq := bsatservice.RmStrDoubles(fileLabel[allloc])
				groupBody = fmt.Sprintf("%v%v\t", groupBody, strings.Join(uniq, " "))
			}

			if len(group) != 0 {
				fileGroup[allloc] = append(fileGroup[allloc], group[gname])
				if len(label) != 0 {
					fileLabel[allloc] = append(fileLabel[allloc], label[gname])
				}
			} else {
				groupBody = "\t"
			}
		}

		// fmt.Printf("Sorting data: Working on %v from %v (%v) \t\t Time:\t%v\n", i, len(*files), gname, t1.Sub(t0))
		// i++
	}

	for _, allloc := range usedLocuses {

		prod := bsatseq.GetProdByName(allloc)

		if bsatstruct.Flag.StatAll == true {
			buffer.WriteString(fmt.Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))

		} else if countNbrOne[allloc] < len(*files) && bsatstruct.Flag.StatAll == false && countNonZeroValues[allloc] > th {
			// fmt.Println(countNbrOne[allloc], allloc)
			// prod := getProductByName(allloc)
			// fmt.Println(strings.Join(filePerLocus[allloc], "\t"))
			buffer.WriteString(fmt.Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))
		}
	}

	// fmt.Println(dndsPerLocus)

	// for _, gname := range genomes {

	// for _, allloc := range usedLocuses {
	// 	prod := getProductByName(allloc)
	// 	// dnds := locInGenome[gname][allloc]

	// 	if countNbrOne[allloc] != len(*files) && *statAll == false {

	// 		buffer.WriteString(fmt.Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))
	// 	} else if *statAll == true {

	// 		buffer.WriteString(fmtdndsPerLocus[allloc] = append(dndsPerLocus[allloc], locInGenome[gname][allloc]).Sprintln(allloc, groupBody, prod, "\t", strings.Join(dndsPerLocus[allloc], "\t"), "\t"))
	// 		// 			buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
	// 	}

	// }

	// }

	headers.WriteString("\n")
	bsatservice.MatrixPrint(headers, buffer, fileOut)

	// headers.WriteString("\n")
	// matrixPrint(headers, buffer, fileOut)

	// fmt.Println(locDNDS)

	// for _, allloc := range allLocuses {

	// 	if countNbrOne[allloc] != len(*files) && *statAll == false {
	// 		buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
	// 	} else if *statAll == true {
	// 		buffer.WriteString(fmt.Sprintln(fmt.Sprintf("%v(%v)\t", allloc, getProductByName(allloc)), strings.Join(locDNDS[allloc], "\t")))
	// 	}

	// }
	// // }

	// headers.WriteString("\n")

	// if buffer.Len() != 0 && headers.Len() != 0 {
	// 	fOut, err := os.Create(fileOut)
	// 	if err != nil {
	// 		log.Fatal("Cannot create file", err)
	// 	}
	// 	defer fOut.Close()
	// 	fmt.Fprintf(fOut, headers.String())
	// 	fmt.Fprintf(fOut, buffer.String())
	// 	fmt.Printf("\n\nWell done!\n")
	// 	// t1 := time.Now()
	// 	// fmt.Printf("Elapsed time: %v", fmtDuration(t1.Sub(t0)))
	// }
	// matrixPrint(headers, buffer, fileOut)
	if bsatstruct.Flag.StatDnDsStat {
		var posGenesHeader strings.Builder
		var posGenesBody strings.Builder
		var posGenesArray = make(map[string][]string)
		uniqPosGenes := bsatservice.RmStrDoubles(positiveGenesList)

		posGenesHeader.WriteString(fmt.Sprintf("genome\t%v\n", strings.Join(uniqPosGenes, "\t")))
		fOut, err := os.Create("dnds_info.txt")

		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer func(fOut *os.File) {
			err := fOut.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOut)

		fOutPos, err := os.Create("dnds_pos_genes.txt")

		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer func(fOutPos *os.File) {
			err := fOutPos.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOutPos)

		_, _ = fmt.Fprint(fOut, fmt.Sprintf("Genome\tNeutral\tNegative\tPositive\tGenesUnderPosSelection\n"))
		for _, fname := range *files {

			_, _ = fmt.Fprint(
				fOut, fmt.Sprintf(
					"%v\t%v\t%v\t%v\t%v\n", fname, countNeutral[fname], countNegative[fname], countPositive[fname], strings.Join(positiveGenes[fname], ",")))

			for _, val := range uniqPosGenes {

				if geneDnDs[fname][val] != "" {
					posGenesArray[fname] = append(posGenesArray[fname], geneDnDs[fname][val])
					// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", fname))
					// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", geneDnDs[fname][val]))
					// } else {
					// 	posGenesArray[fname] = append(posGenesArray[fname], "1.0")
				} else if positiveGenesCheck[val] == 1 {
					posGenesArray[fname] = append(posGenesArray[fname], "NA")
				}
				// posGenesBody.WriteString(fmt.Sprint("\n"))
			}

			// for _, val := range uniqPosGenes {

			// 	// if geneDnDs[fname][val] != 0 {
			// 	// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", fname))
			// 	// 	posGenesBody.WriteString(fmt.Sprintf("%v\t", geneDnDs[fname][val]))
			// 	// }
			// 	// posGenesBody.WriteString(fmt.Sprint("\n"))
			// }

		}

		for _, fname := range *files {
			posGenesBody.WriteString(fmt.Sprintf("%v\t%v\n", fname, strings.Join(posGenesArray[fname], "\t")))
		}

		_, _ = fmt.Fprint(fOutPos, fmt.Sprintf("%v", posGenesHeader.String()))
		_, _ = fmt.Fprint(fOutPos, fmt.Sprintf("%v", posGenesBody.String()))

		fmt.Println("DnDs statistics writen to file dnds_info.txt and dnds_pos_genes.txt")

	}
	// fmt.Println(countNegative, countPositive, countNeutral)
}

func Locus2Matrix(locus string, listOfFiles []string, typeof string) {
	var (
		// coords = regexp.MustCompile(`^(\S+)\W+(\d+)\W+(\d+)`)
		// coordsTwo  = regexp.MustCompile(`^(\w+)\W+(\d+)\W+(\d+)`)
		start, end        int
		altseq, tabMatrix string
		altPostitions     []bsatstruct.AllPositionsInGene
		locmatrix         []string
		locPosLetter      = make(map[int][]string)
	)

	// f, err := os.Open(file) // открываем файл

	// if err != nil {
	// 	fmt.Println(err)

	// }
	// defer f.Close()

	// scanner := bufio.NewScanner(f) //  новый сканер

	// for scanner.Scan() {

	// 	scanTxt := scanner.Text()

	// 	for _, pos := range coords.FindAllStringSubmatch(scanTxt, -1) {
	// 		start, _ = strconv.Atoi(pos[2])
	// 		end, _ = strconv.Atoi(pos[3])
	// fmt.Println(pos)
	// fmt.Println(exGene)
	// exGenes[exGene[0]] = 1
	// exGenes[start] = end
	for _, file := range bsatstruct.ListOfFiles {
		start, end = bsatseq.GetGenePosByName(locus)
		// prod, _ = getProductByPos(start, end)
		altPostitions = bsatvcf.GetAltPos(start, end, file)
		// fmt.Println(altPostitions)
		altseq = bsatseq.MakeAltString(locus, altPostitions)
		locmatrix = strings.Split(altseq, "")
		// fmt.Println(locmatrix)

		// for i := 0; i < len(locmatrix)-1; i++ {
		switch typeof {

		case "nc":
			for i, val := range locmatrix {
				locPosLetter[i] = append(locPosLetter[i], val)

			}
		case "nc_coded":
			for i, val := range locmatrix {
				locPosLetter[i] = append(locPosLetter[i], bsatseq.Nuc2IntCode(val))

			}
		case "binary":
			for i, val := range locmatrix {
				if val == strings.ToUpper(val) {
					locPosLetter[i] = append(locPosLetter[i], "1")
				} else {
					locPosLetter[i] = append(locPosLetter[i], "0")
				}

			}

		}

		// tabMatrix = strings.Join(locmatrix, "\t")
	}

	// result = altStringMatrix{start: start, end: end, locus: locus, tabMatrix: tabMatrix, prod: prod, vcfFile: vcfFile}
	// fmt.Println(len(altPostitions), start, end, locus, prod, altseq)

	// }

	// }

	fmt.Println(strings.Join(listOfFiles, "\t"))

	for i := 0; i < len(locmatrix)-1; i++ {

		tabMatrix = strings.Join(locPosLetter[i], "\t")
		fmt.Println(tabMatrix)
	}

	// return tabMatrix
}

func GetDnDsByLocus(locus string, altPositions []bsatstruct.AllPositionsInGene) (dndsRes []string) {
	var (
		dndsLoc string
	)

	refS := bsatseq.GetGeneSeq(locus)
	altS := bsatseq.MakeAltString(locus, altPositions) // fmt.Println(val, "\n", refS)

	qDnDs := &codon.DnDsQuery{RefSeq: refS, AltSeq: altS, OutChan: make(chan codon.DnDs)}
	go qDnDs.Request()

	dnds := <-qDnDs.OutChan

	close(qDnDs.OutChan)

	if dnds.ND != 0 && dnds.NS != 0 {

		dndsLoc = fmt.Sprintf("%.2f", dnds.DNDS)

	} else {

		dndsLoc = "1.00"

	}
	dndsRes = append(dndsRes, locus)
	dndsRes = append(dndsRes, dndsLoc)

	return dndsRes

}
