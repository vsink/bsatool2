package bsatseq

import (
	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bsatoolMod/src/bsatvcf"
	"bsatoolMod/src/codonPkg"
	"bufio"
	"github.com/pkg/browser"
	"html/template"
	"log"
	"math/rand"
	"net/http"
	"os"
	"regexp"
	"sort"
	"strconv"
	"time"

	"fmt"
	// "log"
	// "math/rand"
	// "os"
	// "sort"
	// "strconv"
	"strings"
	// "time"
)

var (
	tLMN     = "LMN"     // locus:Mutation:NAME
	tPMLN    = "PMLN"    // position:Mutation:locus:NAME
	tPMN     = "PMN"     // position:Mutation:NAME
	tLSAAN   = "LSAAN"   // locus:shortAA:codon:shortAA:name
	tLLAAN   = "LLAAN"   // locus:longAA:codon:longAA:name
	tLCN     = "LCN"     // locus:codon:name
	tSEN     = "SEN"     // start|end:name (any)
	tPNA     = "PNA"     // position~name~antibiotic
	tGLSAAAB = "GLSAAAB" // geneOrLocus_shortAACodonshortAA whiB6_P38L_SM  Rv1258c_V219A_SM
	// tPMNT  = "PMNL"  //position_Ref>Alt{Tab}Name;tag (position, mutation,name, tag )
	// tPOS   = "checkPOS"   // check:apos:name
	// tCODON = "checkCODON" // check:locus:codon:name
	// vcfExt = "*.vcf"
	// pbtmpl      = `{{counters . }}{{ bar . "⎜" "agtc"   "⬤" "TCAG" "⎜"}}{{percent . }}{{rtime . " ETA %s  "}}{{speed .}} `
	lmnRegExp  = `^(\w+)\W+(\d+)(\D)>(\D)\s+(.*)` // LMN Regular expression
	pmlnRegExp = `^(\d+)_(\D)>(\D)\W+(\w+)\W+(.*)`
	pmnRegExp  = `(\d+)_(\D)>(\D)\W+(.*)$`
	// pmntRegExp  = `^(\d+)_(\D)>(\D)(.*)\b(;|:)\b(.*)`
	lsaanRegExp   = `^(\w+)\W+(\D{1})(\d+)(\D{1})\W+(.*)`
	llaanRegExp   = `^(\w+)\W+(\w{3})\W+(\d+)\W+(\w{3})\W+(.*)`
	lcnRegExp     = `^(\w+)\W+codon(\d+)\W+(.*)`
	senRegExp     = `^(\d+)\|(\d+):(\w+)`
	tpnaRegExp    = `^(\d+)~(\S+)~(\w+)`
	glsaaABRegExp = `^(\S+)_(\w)(\d+)(\w)_(\w+)`
)

func ValidateSNP(inString string) bsatstruct.SnpCheckInfo {
	var (
		snpCheck bsatstruct.SnpCheckInfo
	)
	rLMN := regexp.MustCompile(lmnRegExp) // L:1 PiG:2 REF:3 ALT:4 NAME:5
	// LocusMutationName(LMN)
	rPMLN := regexp.MustCompile(pmlnRegExp) // APOS:1 REF:2 ALT:3 L:4 NAME:5
	// PositionMutationLocusName (PMLN)
	rPMN := regexp.MustCompile(pmnRegExp)     // APOS:1 REF:2 ALT:3 NAME:4
	rLSAAN := regexp.MustCompile(lsaanRegExp) // L:1 AA_REF:2 POS:3 AA_ALT:4 NAME:5
	// LocusShortAminoAcidName (LSAAN)
	rLLAAN := regexp.MustCompile(llaanRegExp) // L:1 LAA_REF:2 PiG:3 LAA_ALT:4 NAME:5
	// LocusLongAminoAcidName (LLAAN)
	rLCN := regexp.MustCompile(lcnRegExp)
	// rPMNT := regexp.MustCompile(pmntRegExp) //Pos(1)_Ref(2)>Alt(3){Tab}NAME(4):|;(5)Tag(6)
	rSEN := regexp.MustCompile(senRegExp) // start(1)|end(2):(name)

	rPNA := regexp.MustCompile(tpnaRegExp)        // position~name~antibiotic
	rGLSAAAB := regexp.MustCompile(glsaaABRegExp) // gene(locus)_shortAACodonShortAA

	for _, matchGLSAAAB := range rGLSAAAB.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{Locus: matchGLSAAAB[1], AASref: matchGLSAAAB[2], CodonNbrInG: matchGLSAAAB[3], AASalt: matchGLSAAAB[4], TypeOf: tGLSAAAB,
			AB:  matchGLSAAAB[5],
			Raw: matchGLSAAAB[0]}

	}

	for _, matchLMN := range rLMN.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{Locus: matchLMN[1], PosInGene: matchLMN[2], Ref: matchLMN[3], Alt: matchLMN[4], Name: matchLMN[5], TypeOf: tLMN, Raw: matchLMN[0]}

	}
	for _, matchPMLN := range rPMLN.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{APos: matchPMLN[1], Ref: matchPMLN[2], Alt: matchPMLN[3], Locus: matchPMLN[4], Name: matchPMLN[5], TypeOf: tPMLN, Raw: matchPMLN[0]}

	}
	for _, matchPMN := range rPMN.FindAllStringSubmatch(inString, -1) {

		snpCheck = bsatstruct.SnpCheckInfo{APos: matchPMN[1], Ref: matchPMN[2], Alt: matchPMN[3], Name: matchPMN[4], TypeOf: tPMN, Raw: matchPMN[0]}

	}

	for _, matchLSAAN := range rLSAAN.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{Locus: matchLSAAN[1], AASref: matchLSAAN[2], CodonNbrInG: matchLSAAN[3], AASalt: matchLSAAN[4], Name: matchLSAAN[5], TypeOf: tLSAAN,
			Raw: matchLSAAN[0]}

	}
	for _, matchLLAAN := range rLLAAN.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{Locus: matchLLAAN[1], AALref: matchLLAAN[2], CodonNbrInG: matchLLAAN[3], AALalt: matchLLAAN[4], Name: matchLLAAN[5], TypeOf: tLLAAN,
			Raw: matchLLAAN[0]}

	}
	for _, matchLCN := range rLCN.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{Locus: matchLCN[1], CodonNbrInG: matchLCN[2], Name: matchLCN[3], TypeOf: tLCN, Raw: matchLCN[0]}

	}
	for _, matchSEN := range rSEN.FindAllStringSubmatch(inString, -1) {
		startRange, err1 := strconv.Atoi(matchSEN[1])
		endRange, err2 := strconv.Atoi(matchSEN[2])
		if err1 == nil && err2 == nil {
			snpCheck = bsatstruct.SnpCheckInfo{StartRange: startRange, EndRange: endRange, Tag: matchSEN[3], TypeOf: tSEN}

		} else {
			fmt.Println(err1, err2)
		}

	}
	for _, matchPNA := range rPNA.FindAllStringSubmatch(inString, -1) {
		snpCheck = bsatstruct.SnpCheckInfo{APos: matchPNA[1], Name: matchPNA[2], AB: matchPNA[3], TypeOf: tPNA}

	}

	return snpCheck
}

func Pileup2multifasta(file string, th int, gstrain string, gname string, mkseq bool, minlen int) {
	var (
		nuc                    string
		pos, first, start, end int
		positions              []int
		sequence               []string
		dp                     = 1
		//  foundEnd bool
		fileFormat1 = regexp.MustCompile(`(^\S+)\t(\d+)\W+(\d+)`)
		fileFormat2 = regexp.MustCompile(`(^\S+)\t(\d+)\W+(\w)\W+(\d+)`)

		endFlg                  bool
		resSequence             string
		karyoBuffer, bandBuffer strings.Builder
		fragmentLen             int
	)

	if gstrain == "" {
		gstrain = bsatstruct.GenomeDescription.Strain
	}

	karyoBuffer.WriteString(fmt.Sprintf("chr -  %v caption %v %v %v\n", gstrain, bsatstruct.GenomeDescription.Start, bsatstruct.GenomeDescription.End, gstrain))

	if file != "" {

		file, err := os.Open(file)
		if err != nil {
			log.Fatal(err)
		}
		defer func(file *os.File) {
			err := file.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(file)

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {

			// формат 2 : геном позиция нуклеотид покрытие
			format2Res := fileFormat2.FindAllStringSubmatch(scanner.Text(), -1)

			for _, val2 := range format2Res {
				if val2[4] != "" {

					pos, _ = strconv.Atoi(val2[2])
					dp, _ = strconv.Atoi(val2[4])
					nuc = val2[3]
					if endFlg == true {
						positions = nil
						endFlg = false
					}

					if first == 0 && dp >= th {
						first = pos

					}

					if pos == first+1 && endFlg == false && dp >= th {
						first = pos

						positions = append(positions, pos)
						sequence = append(sequence, nuc)

					} else {

						if len(positions) != 0 {

							start = positions[0] - 1
							end = positions[len(positions)-1]

							switch bsatstruct.Flag.PileupCircos {
							case false:

								if bsatstruct.Flag.PileupShowSeq == true {
									resSequence = strings.Join(sequence, "")

								}
								if start < end && mkseq == false {
									fragmentLen = end - start
									if fragmentLen >= minlen {
										fmt.Printf("%v\t%v\t%v\t%v\t%v\t\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, end-start), resSequence)
									}

								} else if start < end && mkseq == true {
									// here is + added after start
									fragmentLen = end - start
									if fragmentLen >= minlen {
										fmt.Printf(">%v:%v:%v\n%v\n", fmt.Sprintf("%vL%v", gname, end-start), start, end, strings.Join(sequence, ""))
									}

								} else if start > end {
									fragmentLen = start - end
									if fragmentLen >= minlen {
										fmt.Printf("%v\t%v\t%v\t%v\t\n", gstrain, end, start, fmt.Sprintf("%vL%v", gname, start-end))
									}

								}
								endFlg = true
							case true:

								if start < end {

									karyoBuffer.WriteString(
										fmt.Sprintf(
											"band %v %v %v %v %v lblue\n", gstrain, fmt.Sprintf("%vL%v", gname, end-start),
											fmt.Sprintf("%vL%v", gname, end-start), start, end))
									bandBuffer.WriteString(fmt.Sprintf("%v %v %v %v\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, end-start)))
								} else if start > end {
									karyoBuffer.WriteString(
										fmt.Sprintf(
											"band %v %v %v %v %v lblue\n", gstrain, fmt.Sprintf("%vL%v", gname, start-end),
											fmt.Sprintf("%vL%v", gname, start-end), end, start))
									bandBuffer.WriteString(fmt.Sprintf("%v %v %v %v\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, start-end)))
								}
								endFlg = true
								// fmt.Println(buffer.String())

							}
							first = pos
							sequence = nil
						}

					}
				}

			}

			format1Res := fileFormat1.FindAllStringSubmatch(scanner.Text(), -1)

			for _, val := range format1Res {

				pos, _ = strconv.Atoi(val[2])
				dp, _ = strconv.Atoi(val[3])

				if endFlg == true {
					positions = nil
					endFlg = false
				}

				if first == 0 && dp >= th {
					first = pos

				}

				if pos == first+1 && endFlg == false && dp >= th {
					first = pos
					// fmt.Println("yes", strain, end, first, pos)
					positions = append(positions, pos)
				} else {
					if len(positions) != 0 {
						// fmt.Println(positions[0], positions[len(positions)-1], strain)
						start = positions[0] - 1
						end = positions[len(positions)-1]
						if start < end {
							fmt.Printf("%v\t%v\t%v\t%v\t\n", gstrain, start, end, fmt.Sprintf("%vL%v", gname, end-start))
						} else if start > end {
							fmt.Printf("%v\t%v\t%v\t%v\t\n", gstrain, end, start, fmt.Sprintf("%vL%v", gname, start-end))
						}
						endFlg = true

					}
					first = pos

				}

			}

		}

	}

	if karyoBuffer.Len() != 0 && bandBuffer.Len() != 0 {

		var (
			fKaryo = fmt.Sprintf("%v.karyo", file)
			fBand  = fmt.Sprintf("%v.band", file)
		)
		fileKaryo, err := os.Create(fKaryo)
		if err != nil {
			log.Fatal(err)
		}
		defer func(fileKaryo *os.File) {
			err := fileKaryo.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fileKaryo)

		_, _ = fileKaryo.WriteString(karyoBuffer.String())
		fmt.Println(fKaryo, " was successful created!")

		fileBand, err := os.Create(fBand)
		if err != nil {
			log.Fatal(err)
		}
		defer func(fileBand *os.File) {
			err := fileBand.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fileBand)

		_, _ = fileBand.WriteString(bandBuffer.String())
		fmt.Println(fBand, " was successful created!")

	}

}

func ToBED(genes []bsatstruct.GeneInfo, genomeName string) {
	// var (
	// 	genomeName string
	// )

	if genomeName != "" {
		genomeName = bsatstruct.Flag.StatGenomeName
	} else {
		genomeName = bsatstruct.GenomeDescription.Version

	}
	// fmt.Println(genomeName,"!!!")

	for _, g := range genes {
		switch bsatstruct.Flag.StatBedTypeOf {
		case "":
			if bsatstruct.Flag.StatShowAnnotation == true {
				fmt.Printf("%v\t%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus, g.Product)
			} else if bsatstruct.Flag.StatShowAnnotation == false {
				fmt.Printf("%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus)
			}
		case "CDS":
			if g.TypeOf == "CDS" && bsatstruct.Flag.StatShowAnnotation == true {
				fmt.Printf("%v\t%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus, g.Product)
			} else if g.TypeOf == "CDS" && bsatstruct.Flag.StatShowAnnotation == false {
				fmt.Printf("%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus)
			}
		case "IGR":
			if g.TypeOf == "IGR" {
				fmt.Printf("%v\t%v\t%v\t%v\n", genomeName, g.Start, g.End, g.Locus)
			}
		}

	}

}

func ToCircos(genes []bsatstruct.GeneInfo) {
	var (
		// seq        string
		start, end            int
		bandColor, genomeName string

		// gc         float64
		// buffer strings.Builder
	)
	// fmt.Println("---- ideogram.txt ----")
	if bsatstruct.Flag.StatCircosBandColor != "" {
		bandColor = bsatstruct.Flag.StatCircosBandColor
	}
	if bsatstruct.Flag.StatGenomeName != "" {
		genomeName = bsatstruct.Flag.StatGenomeName
	} else {
		genomeName = bsatstruct.GenomeDescription.Strain
	}
	if bsatstruct.Flag.StatCircosTypeOf == "" || bsatstruct.Flag.StatCircosTypeOf == "cds" || bsatstruct.Flag.StatCircosTypeOf == "igr" {
		fmt.Printf(
			"chr -  %v caption %v %v %v\n", genomeName, bsatstruct.GenomeDescription.Start, bsatstruct.GenomeDescription.End,
			bsatstruct.GenomeDescription.Strain)
	}
	for _, g := range genes {
		start = g.Start
		end = g.End
		// seq = getNucFromGenome(start-1, end)
		// fmt.Println(seq)
		if start != 0 && end != 0 {
			// gc, _, _, _ = codon.GcCodonCalc(seq)
			switch bsatstruct.Flag.StatCircosTypeOf {
			case "":
				fmt.Printf("band %v %v %v %v %v %v\n", genomeName, g.Locus, g.Locus, start, end, bandColor)
			case "cds":
				if g.TypeOf == "CDS" {
					fmt.Printf("band %v %v %v %v %v %v\n", genomeName, g.Locus, g.Locus, start, end, bandColor)
				}
			case "gene":
				if g.TypeOf == "CDS" {
					fmt.Printf("%v %v %v %v\n", genomeName, start, end, g.Locus)
				}
			case "all":
				fmt.Printf("%v %v %v %v\n", genomeName, start, end, g.Locus)

			case "igr":
				if g.TypeOf == "IGR" {
					fmt.Printf("band %v %v %v %v %v %v\n", genomeName, g.Locus, g.Locus, start, end, bandColor)
				}
			}

			// buffer.WriteString(fmt.Sprintf("%v %v %v %.2f\n", gInfo.Strain, start, end, gc))
		}
	}

	// fmt.Println("---- histogram.txt ----")
	// fmt.Println(buffer.String())

}

func MakeSeqFasta(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) []bsatstruct.SeqInfo {

	var (
		AllPos, SelectedPos, TempPos []int
		ResSeq                       []bsatstruct.SeqInfo

		uniqueSNP      = make(map[int]int)
		nbrOfSNP       int
		aaAltCoords    = make(map[string]map[int]string)
		aaRefCoords    = make(map[int]string)
		dpMAP          = make(map[int][]int)
		locus          = make(map[int]string)
		prod           = make(map[int]string)
		filesPOS       = make(map[int][]string)
		indel          = make(map[int]int)
		posFN          = make(map[int][]string)
		posFreq        = map[int][]string{}
		snpCount       = make(map[int]int)
		altPercent     int
		altMinMax      []string
		altMin, altMax int
	)

	/*
		длина последовательности
	*/

	if bsatstruct.Flag.AnnSeqLen != 0 {
		nbrOfSNP = bsatstruct.Flag.AnnSeqLen
	}

	/*
		список файлов
	*/

	files := &bsatstruct.ListOfFiles

	for i, file := range *files {
		/*
			создаем map для присваивания aaAltCoords[имя файла][позиция] = val.AltAAShort
		*/

		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan

		snpsChan := make(chan []bsatstruct.SnpInfo)

		go func() {
			snpsChan <- bsatvcf.GetSNPList(file)
		}()
		snpsRes := <-snpsChan

		/*
			создаем map для каждого файла со всеми альтернативными позициями
		*/
		bsatstruct.SnpCacheMap[file] = snpsRes

		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for fname, snps := range bsatstruct.SnpCacheMap {

		for _, val := range snps {
			// DP - глубина прочтения СНИПа
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			// добавляем в массив имена файлов, в соответствии с позицией
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)

			indel[val.APos] = val.Indel

			// if val.DP < bsatstruct.Flag.GbDP {
			// 	// значение 5- пропустить позицию с глубиной меньше указанного значения
			// 	uniqueSNP[val.APos] = 5
			//
			// }

			// exGenes  - гены, которые необходимо исключить из анализа. Загружаются в виде списка командой --exclude-genes
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// значение 2- пропустить данный ген
						uniqueSNP[val.APos] = 2

						continue

						// } else {
						// 	if uniqueSNP[val.APos] != 1 && uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 3 && uniqueSNP[val.APos] != 4 && uniqueSNP[val.
						// 		APos] != 5 && val.Indel == 0 {
						// 		/*
						// 			значение =1 , включаем позицию в формирование последовательности. Главное условие, что данная
						// 			прозиция является не инделом
						// 		*/
						// 		uniqueSNP[val.APos] = 1
						// 		aaAltCoords[fname][val.APos] = val.AltAAShort
						// 		posFN[val.APos] = append(posFN[val.APos], fname)
						//
						// 	} else if val.Indel == 1 {
						// 		/*
						// 			если в позиции находится индел, то пропускаем значение=4
						// 		*/
						// 		uniqueSNP[val.APos] = 4
						//
						// 		continue
						// 	}
					}
				}
			}
			if len(exSNPs) != 0 {

				if exSNPs[val.APos] == 1 {
					/*
						exSNPs - снипы, которые исключаются из анализа. Загружаются командой --exclude-snp
						значение =3 пропускаем позицию
					*/

					uniqueSNP[val.APos] = 3

					continue

					// fmt.Println(val.APos)
				}
			}
			if val.Indel == 1 {
				/*
					если в позиции находится индел, то пропускаем значение=4
				*/
				uniqueSNP[val.APos] = 4

				continue
				// 	}
			}

			if uniqueSNP[val.APos] != 1 && uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 3 && uniqueSNP[val.APos] != 4 {
				uniqueSNP[val.APos] = 1
				aaAltCoords[fname][val.APos] = val.AltAAShort
				aaRefCoords[val.APos] = val.RefAAShort

			}
			if uniqueSNP[val.APos] == 1 {
				if strings.Contains(strings.Join(posFN[val.APos], ""), fname) == false {
					posFN[val.APos] = append(posFN[val.APos], fname)
					snpCount[val.APos] = snpCount[val.APos] + 1
				}

			}

		}

	}

	// fmt.Println(snpCount)
	if bsatstruct.Flag.AnnMaxPos == 0 {
		bsatstruct.Flag.AnnMaxPos = len(*files) - 1
	}
	for key, value := range uniqueSNP {

		if value == 1 && snpCount[key] >= bsatstruct.Flag.AnnMinPos && snpCount[key] <= bsatstruct.Flag.AnnMaxPos {
			altPercent = (snpCount[key] * 100) / len(*files)
			// if len(dpMAP[key]) >= bsatstruct.Flag.GbDP {
			if len(bsatstruct.Flag.AnnAltRange) != 0 {
				altMinMax = strings.Split(bsatstruct.Flag.AnnAltRange, ":")
				if len(altMinMax) != 0 {
					altMin, _ = strconv.Atoi(altMinMax[0])
					altMax, _ = strconv.Atoi(altMinMax[1])
				}
				if altPercent < altMin {
					altPercent = altMin
				}
				if altPercent >= altMin && altPercent <= altMax {
					AllPos = append(AllPos, key)
				}
			} else {
				AllPos = append(AllPos, key)
			}

			// altMinMax = strings.Split(bsatstruct.Flag.AnnAltRange, ":")
			//
			// if len(altMinMax) != 0 {
			// 	altMin, _ = strconv.Atoi(altMinMax[0])
			// 	altMax, _ = strconv.Atoi(altMinMax[1])
			// }
			// if altPercent < altMin {
			// 	altPercent = altMin
			// }
			//
			//
			// if key == 1342096 {
			// 	fmt.Println(posFN[key])
			// }
			if bsatstruct.Flag.GbDebug {
				fmt.Printf("pos:%v count:%v files:%v perсent:%v perc_min:%v perc_max:%v \n", key, snpCount[key], len(*files), altPercent, altMin, altMax)
			}

		}
	}

	sort.Ints(AllPos)

	for _, allpos := range AllPos {
		for _, file := range *files {

			if strings.Contains(strings.Join(posFN[allpos], " "), file) {
				posFreq[allpos] = append(posFreq[allpos], "1")
			} else {
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
		// fmt.Println(allpos, posFreq[allpos], posFN[allpos], strings.Join(posFN[allpos], " "))
	}

	if nbrOfSNP == 0 {
		nbrOfSNP = len(AllPos) - 1
	}

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true && nbrOfSNP != 0 {
		posAlreadyIs := make(map[int]int)
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i

			if posAlreadyIs[rnd] < 1 {
				TempPos = append(TempPos, AllPos[rnd])
				posAlreadyIs[rnd] = posAlreadyIs[rnd] + 1
			}

		}
	} else if randomize == false {

		for i := 0; i <= nbrOfSNP; i++ {

			TempPos = append(TempPos, AllPos[i])
		}
	}

	sort.Ints(TempPos)

	for _, pos := range TempPos {
		var count0, count1 int

		for i := 0; i < len(posFreq[pos]); i++ {
			if posFreq[pos][i] == "0" {
				count0++
			} else if posFreq[pos][i] == "1" {
				count1++
			}
			// fmt.Println(posFreq[pos][i])
		}

		if count0 != 0 {

			// altPercent = (count1 * 100) / len(posFreq[pos])
			// altMinMax = strings.Split(bsatstruct.Flag.AnnAltRange, ":")
			//
			// if len(altMinMax) != 0 {
			// 	altMin, _ = strconv.Atoi(altMinMax[0])
			// 	altMax, _ = strconv.Atoi(altMinMax[1])
			// }
			// if altPercent < altMin {
			// 	altPercent = altMin
			// }
			// if altPercent >= altMin && altPercent <= altMax {
			SelectedPos = append(SelectedPos, pos)
			if bsatstruct.Flag.GbDebug {
				fmt.Printf("pos:%v isRef:%v,isAlt:%v %v legnth_array: %v AltPerc: %v \n", pos, count0, count1, posFreq[pos], len(posFreq[pos]), altPercent)

			}
			// }

		}

	}

	// -pos_file ФЛАГ

	sort.Ints(SelectedPos)

	if bsatstruct.Flag.AnnPosFile != "" && typeof == bsatstruct.NCFlg {
		fOut, err := os.Create(bsatstruct.Flag.AnnPosFile)
		fOutF, err := os.Create(fmt.Sprintf("%v_files.txt", bsatstruct.Flag.AnnPosFile))
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer func(fOut *os.File) {
			err := fOut.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOut)
		defer func(fOutF *os.File) {
			err := fOutF.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOutF)

		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("---NC method----"))
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[pos]apos:indel:dp:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], dpMAP[value], locus[value], prod[value]))
		}
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("---NC method----"))
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[pos]apos:indel:file:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], filesPOS[value], locus[value], prod[value]))
		}

	}

	if ref == true {
		switch typeof {
		case bsatstruct.NCFlg:
			var refBuffer strings.Builder
			refBuffer.WriteString(fmt.Sprintf(">%v\n", "REFERENCE"))
			for _, allpos := range SelectedPos {
				refBuffer.WriteString(bsatservice.GetNucFromGenomePos(allpos))
			}
			ResSeq = append(ResSeq, bsatstruct.SeqInfo{Name: "reference", Seq: refBuffer.String(), UsedPositions: SelectedPos, TypeOfSeq: "NC"})

		}
	}

	for fname, snps := range bsatstruct.SnpCacheMap {

		pos := make(map[int]string)

		var buffer strings.Builder
		// var i int
		var aaSNPpos []int

		buffer.WriteString(fmt.Sprintf(">%v\n", fname))
		switch typeof {
		case bsatstruct.NCFlg:
			for _, val := range snps {
				pos[val.APos] = val.Alt

			}

			for _, allpos := range SelectedPos {

				if pos[allpos] != "" {

					buffer.WriteString(pos[allpos])

				} else {

					buffer.WriteString(bsatservice.GetNucFromGenomePos(allpos))

				}

			}
			ResSeq = append(
				ResSeq, bsatstruct.SeqInfo{Name: fname, Seq: buffer.String(), UsedPositions: SelectedPos, TypeOfSeq: "AA"})

		case bsatstruct.AAFlg:

			for _, allpos := range SelectedPos {
				if aaAltCoords[fname][allpos] == "" {
					buffer.WriteString(strings.ToLower(aaRefCoords[allpos]))
					aaSNPpos = append(aaSNPpos, allpos)
				} else if aaAltCoords[fname][allpos] != "" {
					aaSNPpos = append(aaSNPpos, allpos)
					buffer.WriteString(aaAltCoords[fname][allpos])
				}

			}

			ResSeq = append(ResSeq, bsatstruct.SeqInfo{Name: fname, Seq: buffer.String(), UsedPositions: aaSNPpos})

		}
	}

	return ResSeq
}

func MakeAlign(verbose bool) {

	// type posPerFile struct {
	// 	Fname string
	// 	Pos   []map[int]string
	// }

	var (
		uniqueSNP = make(map[int]int)

		aaAltCoords = make(map[string]map[int]string)

		dpMAP    = make(map[int][]int)
		locus    = make(map[int]string)
		prod     = make(map[int]string)
		filesPOS = make(map[int][]string)
		indel    = make(map[int]int)

		lineSeq []string
	)

	files := &bsatstruct.ListOfFiles

	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		// qSNP := bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery), Print: verbose}
		// go qSNP.Request()
		// // snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		// snpRes := <-qSNP.OutChan
		// bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo

		snpsChan := make(chan []bsatstruct.SnpInfo)

		go func() {
			snpsChan <- bsatvcf.GetSNPList(file)
		}()
		snpsRes := <-snpsChan
		bsatstruct.SnpCacheMap[file] = snpsRes

		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for fname, snps := range bsatstruct.SnpCacheMap {
		fmt.Printf("\n>%v\n", fname)
		snpPos := make(map[int]string)
		for _, val := range snps {
			// DP - глубина прочтения СНИПа
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)

			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)
			indel[val.APos] = val.Indel
			if val.DP < bsatstruct.Flag.GbDP {
				uniqueSNP[val.APos] = 2
			}
			if val.Indel == 0 {

				snpPos[val.APos] = val.Alt

			} else {

				snpPos[val.APos] = val.IndelAlt
			}

		}

		for i := bsatstruct.GenomeDescription.Start; i < bsatstruct.GenomeDescription.End; i++ {

			// lenCount++
			// fmt.Println(i,getNucFromGenomePos(i),snpPos[i])
			if len(snpPos[i]) != 0 {
				// lineSeq.WriteString(snpPos[i])
				lineSeq = append(lineSeq, snpPos[i])
				// fmt.Println(snpPos[i])
			} else {
				lineSeq = append(lineSeq, bsatservice.GetNucFromGenomePos(i))
				// fmt.Println(getNucFromGenomePos(i))
			}
			if len(lineSeq) == 70 {
				// seq.WriteString(lineSeq.String())
				// fmt.Println(lineSeq.String())
				seq := strings.Join(lineSeq, "")
				fmt.Printf("%v\n", seq)
				lineSeq = nil
				// seq.WriteString("\n")

			}
		}
	}

	// }()
	// wg.Wait()

	// fmt.Println(posPerFile)
}

func MakeSeqBin(typeof string, verbose bool, exGenes map[int]int, exSNPs map[int]int, randomize bool, dp int) {

	var (
		AllPos, SelectedPos []int

		// passSNP = make(map[string]int)
		uniqSNP     = make(map[int]int)
		nbrOfSNP    int
		aaAltCoords = make(map[string]map[int]string)
		aaRefCoords = make(map[int]string)
		dpMAP       = make(map[int][]int)
		locus       = make(map[int]string)
		prod        = make(map[int]string)
		indel       = make(map[int]int)
		filesPOS    = make(map[int][]string)
		nexusTaxa   = make(map[string][]string)
		nChar       int
		nTax        int
	)

	if bsatstruct.Flag.AnnSeqLen != 0 {
		nbrOfSNP = bsatstruct.Flag.AnnSeqLen
	}

	files := &bsatstruct.ListOfFiles
	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		// qSNP := bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery), Print: verbose}
		// go qSNP.Request()
		// snpRes := <-qSNP.OutChan
		// bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
		snpsChan := make(chan []bsatstruct.SnpInfo)

		go func() {
			snpsChan <- bsatvcf.GetSNPList(file)
		}()
		snpsRes := <-snpsChan
		bsatstruct.SnpCacheMap[file] = snpsRes
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for fname, snps := range bsatstruct.SnpCacheMap {

		// fmt.Println(tmp)
		for _, val := range snps {
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			indel[val.APos] = val.Indel
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)

			if val.DP < bsatstruct.Flag.GbDP {
				uniqSNP[val.APos] = 2
			}
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {

						uniqSNP[val.APos] = 2

						continue
					} else if exSNPs[val.APos] == 1 {

						uniqSNP[val.APos] = 2

						continue

					} else {
						if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 && val.DP >= dp {
							uniqSNP[val.APos] = 1
							aaAltCoords[fname][val.APos] = val.AltAAShort

						}
					}
				}
			} else {
				if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
					uniqSNP[val.APos] = 1
					aaAltCoords[fname][val.APos] = val.AltAAShort
					aaRefCoords[val.APos] = val.RefAAShort

				}

			}

		}
	}

	for key, value := range uniqSNP {

		if value == 1 {

			if len(dpMAP[key]) >= bsatstruct.Flag.AnnMinPos {
				AllPos = append(AllPos, key)

			}

		}

	}

	sort.Ints(AllPos)

	if nbrOfSNP == 0 {
		nbrOfSNP = len(AllPos) - 1
	}

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true && bsatstruct.Flag.AnnSeqLen != 0 {
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
			if len(AllPos) != 0 {
				SelectedPos = append(SelectedPos, AllPos[i])
			}

		}
	}

	sort.Ints(SelectedPos)

	if bsatstruct.Flag.AnnPosFile != "" {
		fOut, err := os.Create(bsatstruct.Flag.AnnPosFile)
		fOutF, err := os.Create(fmt.Sprintf("%v_files.txt", bsatstruct.Flag.AnnPosFile))
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer func(fOut *os.File) {
			err := fOut.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOut)
		defer func(fOutF *os.File) {
			err := fOutF.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOutF)
		// var posArr []string
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("---Binary method----"))
		_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[pos]apos:indel:dp:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOut, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], dpMAP[value], locus[value], prod[value]))
		}
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("---Binary method----"))
		_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[pos]apos:indel:file:loc:prod"))
		for i, value := range SelectedPos {

			_, _ = fmt.Fprintln(fOutF, fmt.Sprintf("[%v]%v:%v:%v:%v:%v\t", i, value, indel[value], filesPOS[value], locus[value], prod[value]))
		}

	}
	nTax = len(bsatstruct.SnpCacheMap)
	nChar = len(SelectedPos)
	// fmt.Println(AllPos)
	if bsatstruct.Flag.AnnOutFormat == "phylip" {
		fmt.Printf("%v\t%v\n", nTax, nChar)
	}
	for fname, snps := range bsatstruct.SnpCacheMap {

		pos := make(map[int]string)

		var buffer strings.Builder

		switch typeof {
		case "Binary":

			if bsatstruct.Flag.AnnOutFormat == "phylip" {
				for _, val := range snps {
					pos[val.APos] = val.Alt

				}
				for _, allpos := range SelectedPos {
					// posCount[allpos] = posCount[allpos] + 1
					if pos[allpos] != "" {

						buffer.WriteString("1")
					} else {

						buffer.WriteString("0")
					}

				}

				fmt.Println(fname, "\t", buffer.String())

			} else if bsatstruct.Flag.AnnOutFormat == "nexus" {

				for _, val := range snps {
					pos[val.APos] = val.Alt

				}
				for _, allpos := range SelectedPos {
					// posCount[allpos] = posCount[allpos] + 1
					if pos[allpos] != "" {

						nexusTaxa[fname] = append(nexusTaxa[fname], "1")
					} else {

						nexusTaxa[fname] = append(nexusTaxa[fname], "0")
					}

				}

			}

		}
	}

	if len(nexusTaxa) != 0 {
		var taxNbr int
		// fmt.Println(nexusTaxa)
		fmt.Printf("#nexus \n\nBEGIN Taxa;\nDIMENSIONS\nntax=%v;\nTAXLABELS\n", nTax)
		for key := range nexusTaxa {
			taxNbr++
			fmt.Printf("[%v]\t'%v'\n", taxNbr, key)
		}
		fmt.Printf(
			";\nEND; [Taxa]\n\nBEGIN Characters;\nDIMENSIONS nchar=%v;\nFORMAT\n\tdatatype=STANDARD\n\tmissing=?\n\tgap=-\n\tsymbols=\"01\"\n\t labels=left\n\ttranspose=no\n\tinterleave=no\n;\nMATRIX\n",
			nChar)
		for key, val := range nexusTaxa {

			fmt.Printf("'%v'\t%v\n", key, strings.Join(val, ""))
		}

		fmt.Printf("\n;\nEnd;\n")
	}

}

func MakeSeqNex(typeof string, verbose bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) {

	var (
		AllPos, SelectedPos, TempPos []int
		// passSNP = make(map[string]int)
		uniqueSNP      = make(map[int]int)
		nbrOfSNP       int
		aaAltCoords    = make(map[string]map[int]string)
		aaRefCoords    = make(map[int]string)
		dpMAP          = make(map[int][]int)
		locus          = make(map[int]string)
		prod           = make(map[int]string)
		filesPOS       = make(map[int][]string)
		indel          = make(map[int]int)
		posFN          = make(map[int][]string)
		posFreq        = map[int][]string{}
		altPercent     int
		altMinMax      []string
		altMin, altMax int
		nexusTaxa      = make(map[string][]string)
		nChar          int
		nTax           int
	)

	if bsatstruct.Flag.AnnSeqLen != 0 {
		nbrOfSNP = bsatstruct.Flag.AnnSeqLen
	}

	files := &bsatstruct.ListOfFiles
	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		// qSNP := bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery), Print: verbose}
		// go qSNP.Request()
		// // snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		// snpRes := <-qSNP.OutChan
		// bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
		snpsChan := make(chan []bsatstruct.SnpInfo)

		go func() {
			snpsChan <- bsatvcf.GetSNPList(file)
		}()
		snpsRes := <-snpsChan
		bsatstruct.SnpCacheMap[file] = snpsRes
		if verbose == true {
			fmt.Printf("Reading files:  %v from %v \r", i+1, len(*files))
		}

	}

	for fname, snps := range bsatstruct.SnpCacheMap {

		for _, val := range snps {
			// DP - глубина прочтения СНИПа
			dpMAP[val.APos] = append(dpMAP[val.APos], val.DP)
			// indel[val.APos] = val.Indel
			locus[val.APos] = val.Locus
			prod[val.APos] = val.Product
			filesPOS[val.APos] = append(filesPOS[val.APos], fname)
			indel[val.APos] = val.Indel
			if val.DP < bsatstruct.Flag.GbDP {
				uniqueSNP[val.APos] = 2
			}

			// exGenes  - гены, которые необходимо исключить из анализа. Загружаются в виде списка командой --exclude-genes
			if len(exGenes) != 0 {

				for key, value := range exGenes {
					if val.APos >= key && val.APos <= value {
						// значение 2- пропустить данный ген
						uniqueSNP[val.APos] = 2

						continue
					} else if exSNPs[val.APos] == 1 {
						/*
							exSNPs - снипы, которые исключаются из анализа. Загружаются командой --exclude-snp
							значение =2 пропускаем позицию
						*/
						uniqueSNP[val.APos] = 2
						continue

					} else {
						if uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 1 && val.Indel == 0 {
							/*
								значение =1 , включаем позицию в формирование последовательности. Главное условие, что данная
								прозиция является не инделом
							*/
							uniqueSNP[val.APos] = 1
							aaAltCoords[fname][val.APos] = val.AltAAShort
							posFN[val.APos] = append(posFN[val.APos], fname)

						} else if val.Indel == 1 {
							/*
								если в позиции находится индел, то пропускаем
							*/
							uniqueSNP[val.APos] = 2
							continue
						}
					}
				}
			} else {
				if uniqueSNP[val.APos] != 2 && uniqueSNP[val.APos] != 1 && val.Indel == 0 {
					uniqueSNP[val.APos] = 1
					aaAltCoords[fname][val.APos] = val.AltAAShort
					aaRefCoords[val.APos] = val.RefAAShort
				}
				if uniqueSNP[val.APos] == 1 {
					posFN[val.APos] = append(posFN[val.APos], fname)
				}

			}

		}
	}

	for key, value := range uniqueSNP {

		if value == 1 {
			if len(dpMAP[key]) >= bsatstruct.Flag.AnnMinPos {
				AllPos = append(AllPos, key)

			}

		}

	}

	sort.Ints(AllPos)
	for _, allpos := range AllPos {

		for _, file := range *files {
			if strings.Contains(strings.Join(posFN[allpos], " "), file) {
				posFreq[allpos] = append(posFreq[allpos], "1")

			} else {
				// fmt.Println(allpos, 0, file)
				posFreq[allpos] = append(posFreq[allpos], "0")
			}

		}
	}

	if nbrOfSNP == 0 {
		nbrOfSNP = len(AllPos) - 1
	}

	if nbrOfSNP > len(AllPos) {
		nbrOfSNP = len(AllPos) - 1
	}

	if randomize == true && nbrOfSNP != 0 {
		posAlreadyIs := make(map[int]int)
		rand.Seed(time.Now().UnixNano())
		for i := 1; i <= nbrOfSNP; i++ {
			rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])
			if posAlreadyIs[rnd] < 1 {
				TempPos = append(TempPos, AllPos[rnd])
				posAlreadyIs[rnd] = posAlreadyIs[rnd] + 1
			}

		}
	} else if randomize == false {

		for i := 0; i <= nbrOfSNP; i++ {

			TempPos = append(TempPos, AllPos[i])
		}
	}

	sort.Ints(TempPos)

	for _, pos := range TempPos {
		var count0, count1 int
		for i := 0; i < len(posFreq[pos]); i++ {
			if posFreq[pos][i] == "0" {
				count0++
			} else if posFreq[pos][i] == "1" {
				count1++
			}
		}
		if count0 != 0 {

			altPercent = (count1 * 100) / len(posFreq[pos])
			altMinMax = strings.Split(bsatstruct.Flag.AnnAltRange, ":")
			if len(altMinMax) != 0 {
				altMin, _ = strconv.Atoi(altMinMax[0])
				altMax, _ = strconv.Atoi(altMinMax[1])
			}

			if altPercent >= altMin && altPercent <= altMax {
				SelectedPos = append(SelectedPos, pos)
				if bsatstruct.Flag.GbDebug {
					fmt.Printf("pos:%v isRef:%v,isAlt:%v %v legnth_array: %v AltPerc: %v \n", pos, count0, count1, posFreq[pos], len(posFreq[pos]), altPercent)

				}
			}

		}

	}

	// -pos_file ФЛАГ

	sort.Ints(SelectedPos)
	for fname, snps := range bsatstruct.SnpCacheMap {

		pos := make(map[int]string)
		var buffer strings.Builder

		buffer.WriteString(fmt.Sprintf(">%v\n", fname))
		switch typeof {
		case bsatstruct.NCFlg:
			for _, val := range snps {
				pos[val.APos] = val.Alt

			}
			nTax = len(bsatstruct.SnpCacheMap)
			nChar = len(SelectedPos)
			for _, allpos := range SelectedPos {

				if pos[allpos] != "" {

					nexusTaxa[fname] = append(nexusTaxa[fname], pos[allpos])

				} else {

					nexusTaxa[fname] = append(nexusTaxa[fname], bsatservice.GetNucFromGenomePos(allpos))
				}

			}

		}

	}
	if len(nexusTaxa) != 0 {

		fmt.Printf("#NEXUS\nBegin data;\nDimensions ntax=%v nchar=%v\nFormat datatype=dna missing=N gap=-;\n", nTax, nChar)

		fmt.Printf("matrix\n")
		for key, val := range nexusTaxa {

			fmt.Printf("%v  \n%v\n", key, strings.Join(val, ""))
		}
		fmt.Printf(";\nEnd;\n")
	}

}

func CheckVCFFormat(file string) bool {

	var (
		lineCount int
		res       bool
	)

	f, err := os.Open(file) // открываем файл

	if err != nil {
		log.Fatal(file, err.Error())

	}
	defer func(f *os.File) {
		err := f.Close()
		if err != nil {
			log.Fatal(file, err.Error())
		}
	}(f)

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {
		lineCount++
		if strings.Contains(scanner.Text(), "##fileformat=VCFv") && lineCount == 1 {
			res = true
			break
		}

	}
	return res
}

func ParserBulkVCF(withFilenames bool) {

	files := &bsatstruct.ListOfFiles
	for i, file := range *files {

		snpsChan := make(chan []bsatstruct.SnpInfo)

		go func() {
			snpsChan <- bsatvcf.GetSNPList(file)
		}()
		snps := <-snpsChan

		if bsatstruct.Flag.GbVerbose == true {
			fmt.Printf(" File %v (%v from %v) was succefull annotated. \r", file, i+1, len(*files))
		}
		if withFilenames == true {
			fmt.Printf("\n\n%v:\n\n", file)

		}
		// if *annShowFileName == false {
		for _, val := range snps {
			val.FName = file
			bsatservice.PrintInConsole(val, bsatstruct.Flag.GbIGR)
		}

	}

}

func Pos2Seq(coordFile string, targetFile string) bsatstruct.AltStringResult {
	var (
		coords = regexp.MustCompile(`^(\S+)\W+(\d+)\W+(\d+)`)

		start, end          int
		locus, prod, altseq string
		altPostitions       []bsatstruct.AllPositionsInGene
		result              bsatstruct.AltStringResult
	)

	f, err := os.Open(coordFile) // открываем файл

	if err != nil {
		fmt.Println(err)

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

		for _, pos := range coords.FindAllStringSubmatch(scanTxt, -1) {
			start, _ = strconv.Atoi(pos[2])
			end, _ = strconv.Atoi(pos[3])

			locus, _ = bsatservice.GetGeneNameByPos(start, end)
			prod, _ = bsatservice.GetProdByPos(start, end)
			altPostitions = GetAltPosFromFile(start, end, targetFile)

			altseq = MakeAltStringByPos(start, end, altPostitions)

			result = bsatstruct.AltStringResult{Start: start, End: end, Locus: locus, AltSeq: altseq, Prod: prod, VcfFile: targetFile}

		}

	}
	return result
}

func NucWebServer(port string, exGenes map[int]int, exSNPs map[int]int) {
	/*

	 */

	seq := MakeSeqFasta(bsatstruct.NCFlg, bsatstruct.Flag.GbVerbose, bsatstruct.Flag.AnnMakeSeqRef, exGenes, exSNPs, bsatstruct.Flag.GbRandomize)
	var htmlTemplate = `
<!DOCTYPE html>
<html>
<head>
<style>
.col {
word-wrap: break-word; /* Перенос слов */
}
</style>
</head>
<body>
<div class="col">
{{range $element := .}}

<p>{{.Seq}}</p>
		{{end}}
</div>
</body>


`

	t := template.New("t")
	t, err := t.Parse(htmlTemplate)
	if err != nil {
		panic(err)
	}

	http.HandleFunc(
		"/", func(w http.ResponseWriter, r *http.Request) {

			err = t.Execute(w, seq)
			if err != nil {
				panic(err)
			}

			go func() {
				defer os.Exit(0)
			}()
		})

	_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

	locPort := fmt.Sprintf(":%v", port)

	_ = http.ListenAndServe(locPort, nil)

}

func GetAllPositions(exGenes map[int]int, exSNPs map[int]int) (allPos []int) {
	var (
		AllPosUnsort []int
		uniqSNP      = make(map[int]int)
		passedSNP    int
	)
	if len(bsatstruct.SnpCacheMap) != 0 {
		for _, snps := range bsatstruct.SnpCacheMap {

			for _, val := range snps {

				if len(exGenes) != 0 {

					for key, value := range exGenes {
						if val.APos >= key && val.APos <= value {

							uniqSNP[val.APos] = 2 // 2-EXCLUDED

							continue
						} else if exSNPs[val.APos] == 1 {

							uniqSNP[val.APos] = 2 // 2-EXCLUDED

							continue

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

				}

				if uniqSNP[val.APos] == 1 {
					AllPosUnsort = append(AllPosUnsort, val.APos)
				} else {
					passedSNP = passedSNP + 1
				}
			}

		}
		allPos = bsatservice.GenUnique(AllPosUnsort)

		sort.Ints(allPos)

	}
	if bsatstruct.Flag.GbDebug {
		fmt.Println("passed : ", passedSNP, " snps")
	}
	return allPos
}

func GetAltPosFromFile(start int, end int, vcfFile string) []bsatstruct.AllPositionsInGene {
	var (
		snps         []bsatstruct.SnpInfo
		altPositions []bsatstruct.AllPositionsInGene
	)

	snpsChan := make(chan []bsatstruct.SnpInfo)

	go func() {
		snpsChan <- bsatvcf.GetSNPList(vcfFile)
	}()
	snps = <-snpsChan

	for _, val := range snps {

		if start >= val.Start && end <= val.End && val.TypeOf == "CDS" {
			altPositions = append(altPositions, bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})

		}
	}

	return altPositions
}

func Snp2SeqByLoc(locus string, vcfFile string) bsatstruct.AltStringResult {
	var (
		start, end    int
		prod, altseq  string
		altPostitions []bsatstruct.AllPositionsInGene
		result        bsatstruct.AltStringResult
	)

	start, end = bsatservice.GetGenePosByName(locus)
	prod, _ = bsatservice.GetProdByPos(start, end)
	altPostitions = GetAltPosFromFile(start, end, vcfFile)

	altseq = MakeAltString(locus, altPostitions, false)

	result = bsatstruct.AltStringResult{Start: start, End: end, Locus: locus, AltSeq: altseq, Prod: prod, VcfFile: vcfFile}

	return result
}

func GetSeqRange(PosRange string, noseq bool) bsatstruct.RangePosInfo {

	var (
		rangeParser = regexp.MustCompile(`\b(\d+)\b\W+\b(\d+)\b`)

		result bsatstruct.RangePosInfo
	)
	PosRange = strings.Replace(PosRange, ".", "", -1)

	for _, val := range rangeParser.FindAllStringSubmatch(PosRange, -1) {

		if val[1] != "" && val[2] != "" {
			startR, _ := strconv.Atoi(val[1])
			endR, _ := strconv.Atoi(val[2])
			if startR < endR {
				bsatstruct.ResStart = startR
				bsatstruct.ResEnd = endR
			} else {
				bsatstruct.ResStart = endR
				bsatstruct.ResEnd = startR
			}

		}

	}

	if bsatstruct.ResStart != 0 && bsatstruct.ResEnd != 0 {

		if noseq == true {
			gname, _ := bsatservice.GetGeneNameByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			prod, _ := bsatservice.GetProdByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			result = bsatstruct.RangePosInfo{Start: bsatstruct.ResStart + 1, End: bsatstruct.ResEnd, GeneName: gname, Prod: prod}

		} else if noseq == false {
			seq := bsatservice.GetNucFromGenome(bsatstruct.ResStart-1, bsatstruct.ResEnd)
			gname, _ := bsatservice.GetGeneNameByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			prod, _ := bsatservice.GetProdByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			result = bsatstruct.RangePosInfo{Start: bsatstruct.ResStart + 1, End: bsatstruct.ResEnd, GeneName: gname, Len: len(seq), Seq: seq, Prod: prod}

		}

	} else {
		fmt.Println("Range positions is not valid")
	}

	return result

}

func MakeAltString(locus string, positions []bsatstruct.AllPositionsInGene, nucAsDots bool) string {
	// var lStart, lEnd int
	var (
		seqSplit []string
		seq      string
	)

	seq = GetGeneSeq(locus)

	for _, nuc := range seq {
		if nucAsDots == false {
			seqSplit = append(seqSplit, string(nuc))
		} else {
			seqSplit = append(seqSplit, ".")
		}

	}

	for _, val := range positions {
		// fmt.Println(val.pos, val.ref, ">", val.alt, seqSplit[val.pos-1], ">", val.alt)
		if len(seqSplit) != 0 {
			seqSplit[val.Pos-1] = strings.ToUpper(val.Alt)

		}

	}

	return strings.Join(seqSplit, "")

}

func MakeAltStringByPos(start int, end int, positions []bsatstruct.AllPositionsInGene) string {
	// var lStart, lEnd int
	var (
		seqSplit []string
		seq      string
	)

	seq = bsatservice.GetNucFromGenome(start, end)

	for _, nuc := range seq {
		seqSplit = append(seqSplit, string(nuc))

	}

	for _, val := range positions {

		if len(seqSplit) != 0 {
			seqSplit[val.Pos-1] = val.Alt
		}

	}

	return strings.Join(seqSplit, "")

}

func GetGeneSeq(locus string) string {
	var seq string
	start, end := bsatservice.GetGenePosByName(locus)
	if start != 0 && end != 0 {
		dir := bsatservice.GetDirectByName(locus)
		if dir == "f" {
			seq = bsatservice.GetNucFromGenome(start-1, end)
		} else if dir == "r" {
			seq = bsatservice.GetRevComplement(bsatservice.GetNucFromGenome(start-1, end))

		}
	} else {
		fmt.Println("Locus not found [getGeneSequence func]", locus)
	}
	return seq
}

func GetInfoByLocus() {
	var (
		seq, direction              string
		asFasta                     bool
		flankSeqLeft, flankSeqRight string
		flankInfo                   strings.Builder
		colorRed                    = "\x1b[31;1m"
	)
	start, end := bsatservice.GetGenePosByName(bsatstruct.Flag.InfoLocus)
	// fmt.Println(start, end)
	if bsatstruct.Flag.InfoLocus != "" {
		if start != 0 && end != 0 {
			switch bsatstruct.Flag.InfoShowAs {
			case "", "direct":
				seq = bsatservice.GetNucFromGenome(start-1, end)
				direction = "in forward direction"
			case "gene":
				seq = GetGeneSeq(bsatstruct.Flag.InfoLocus)
				direction = "in gene direction"
			case "fasta":
				asFasta = true
				if bsatstruct.Flag.InfoFlankLeft != 0 {
					// flankSeqLeft = bsatservice.GetNucFromGenome((start-bsatstruct.Flag.InfoFlankLeft)-1, end)
					start = start - bsatstruct.Flag.InfoFlankLeft
					flankInfo.WriteString(fmt.Sprintf("+%v_left ", bsatstruct.Flag.InfoFlankLeft))
					// fmt.Println("*", flankSeqLeft, "LEFT")
				}
				if bsatstruct.Flag.InfoFlankRight != 0 {
					flankSeqRight = bsatservice.GetNucFromGenome(end, end+bsatstruct.Flag.InfoFlankRight)
					// fmt.Println("*", flankSeqRight, "RIGHT")
					flankInfo.WriteString(fmt.Sprintf("+%v_right ", bsatstruct.Flag.InfoFlankRight))
				}

				seq = bsatservice.GetNucFromGenome(start-1, end)
				if bsatstruct.Flag.InfoMask != "" {
					posRefAlt := strings.Split(bsatstruct.Flag.InfoMask, ":")
					locPosInGene, _ := strconv.Atoi(posRefAlt[0])
					if locPosInGene != 0 && locPosInGene < len(seq) {
						splitSeq := strings.Split(seq, "")
						if strings.ToUpper(splitSeq[locPosInGene-1]) == strings.ToUpper(posRefAlt[1]) {
							if bsatstruct.Flag.InfoIUPAc == false {

								splitSeq[locPosInGene-1] = fmt.Sprintf("[%v/%v]", posRefAlt[1], posRefAlt[2])
							} else {
								locUIPAC := bsatservice.GetIUPAC(fmt.Sprintf("%v%v", posRefAlt[1], posRefAlt[2]))
								splitSeq[locPosInGene-1] = fmt.Sprintf("\n%v\n", locUIPAC)
							}

						} else {
							splitSeq[locPosInGene-1] = fmt.Sprintf("[%v/%v]", posRefAlt[1], posRefAlt[2])
							fmt.Println("Masked reference nucleotide does not match to reference!")
						}

						seq = strings.Join(splitSeq, "")

						if bsatstruct.Flag.InfoFlankRight != 0 || bsatstruct.Flag.InfoFlankLeft != 0 {
							seq = fmt.Sprintf("%v%v%v", flankSeqLeft, seq, flankSeqRight)
						}

					}
				}

				direction = "in forward direction"
			}

			prod := bsatservice.GetProdByName(bsatstruct.Flag.InfoLocus)
			note := bsatservice.GetNoteByName(bsatstruct.Flag.InfoLocus)
			goa := bsatservice.GetGOAByName(bsatstruct.Flag.InfoLocus)

			// вывод последовательности кодонов заданной ключом --codons=Start:End
			if bsatstruct.Flag.InfoCodons != "" {

				seqByCodons := make([]string, (len(seq)/3)+1)
				codonNbr := 1
				var codonBuffer strings.Builder
				for _, c := range seq {
					codonBuffer.WriteString(string(c))
					if codonBuffer.Len() == 3 {

						seqByCodons[codonNbr] = codonBuffer.String()

						codonBuffer.Reset()
						codonNbr = codonNbr + 1
					}
				}
				codonsCoords := strings.Split(bsatstruct.Flag.InfoCodons, ":")
				startCodon, _ := strconv.Atoi(codonsCoords[0])
				endCodon, _ := strconv.Atoi(codonsCoords[1])
				for i := startCodon; i <= endCodon; i++ {
					if bsatstruct.Flag.GbVerbose == false {
						codonBuffer.WriteString(seqByCodons[i])
					} else {

						codonBuffer.WriteString(fmt.Sprintf("[%v]%v", i, seqByCodons[i]))
					}

				}

				fmt.Printf(">%v %v [%v:%v codons %v]\n%v\n", bsatstruct.Flag.InfoLocus, prod, startCodon, endCodon, direction, codonBuffer.String())

			} else {

				if asFasta == false {
					gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)
					fmt.Printf(
						">%v %v (%v-%v %v)\n%v\n-----------------\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n-----------------\n%v\ngoa:%v\n",
						bsatstruct.Flag.InfoLocus,
						prod,
						start, end, direction, seq, gc, gc1, gc2, gc3, len(seq), note, goa)
				} else {

					fmt.Printf(">%v %v (%v-%v %v) %v \n%v\ngoa:%v\n", bsatstruct.Flag.InfoLocus, prod, start, end, direction, flankInfo.String(), seq, goa)
				}

			}
		} else {
			fmt.Println(colorRed, bsatstruct.Flag.InfoLocus, " not found!")
		}
	} else if bsatstruct.Flag.InfoRanges != "" {
		coords := strings.Split(bsatstruct.Flag.InfoRanges, ":")

		start, _ := strconv.Atoi(coords[0])
		end, _ := strconv.Atoi(coords[1])
		locname, _ := bsatservice.GetGeneNameByPos(start, end)
		if start != 0 && end != 0 {
			seq := bsatservice.GetNucFromGenome(start-1, end)
			// gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)
			// prod := getProductByName(*infoLocus)
			// fmt.Printf("%v (%v-%v)\n%v\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n", prod, start, end, seq, gc, gc1, gc2, gc3, len(seq))
			// fmt.Println(seq, gc, gc1, gc2, gc3, prod)
			fmt.Printf(">%v|%v|%v (%v:%v)\n%v\n", bsatstruct.GenomeDescription.Strain, bsatstruct.GenomeDescription.Version, locname, start, end, seq)

		}

	}
}

func TestGeneInfo(genes []bsatstruct.GeneInfo) {
	var (
		prod              string
		start, end        int
		gc, gc1, gc2, gc3 float64
	)
	for _, g := range genes {

		start = g.Start
		end = g.End

		if start != 0 && end != 0 && g.TypeOf == "CDS" {
			seq := GetGeneSeq(g.Locus)
			prod = bsatservice.GetProdByName(g.Locus)
			gc, gc1, gc2, gc3 = codon.GcCodonCalc(seq)
			fmt.Printf("%v:%v\t%v\t%v\t%.2f\t%.2f\t%.2f\t%.2f\n", start, end, g.Locus, prod, gc, gc1, gc2, gc3)

		}
	}
}

// func MakeChangedSeq(startPos int, endPos int, posToChange int, nucToChange string, flankLeft int, flankRight int) {
// 	var (
// 		// seq, direction              string
// 		seq string
// 		// asFasta                     bool
// 		flankSeqLeft, flankSeqRight string
// 		// flankInfo                   strings.Builder
// 		// colorRed                    = "\x1b[31;1m"
// 		refNuc string
// 	)
//
// 	if startPos != 0 && endPos != 0 && posToChange != 0 {
// 		// dir := bsatservice.GetDirectByName(locus)
//
// 		/*
// 				flankLeft........StartPos......PosToChange.....endPos........flankRight
// 				-------------->------------------->------------------>------------->
// 				seq[startPos-flankLeft]........seq[ALtNuc].......seq[endPos+flankRight]
// 			start-------------------------------------------------------------------> end
// 		*/
//
// 		if startPos < endPos {
// 			fmt.Println("forward")
// 			seq = bsatservice.GetNucFromGenome(startPos-1, endPos)
// 			refNuc = bsatservice.GetNucFromGenomePos(posToChange)
// 			if flankLeft != 0 {
// 				flankSeqLeft = bsatservice.GetNucFromGenome((startPos-flankLeft)-1, startPos)
// 			}
// 			if flankRight != 0 {
// 				flankSeqRight = bsatservice.GetNucFromGenome(endPos-1, endPos+flankRight)
// 			}
// 		} else if startPos > endPos {
// 			fmt.Println("reverse")
// 			/*
// 				flankRight.......StartPos.......PosToChange.....endPos....flankLeft
// 				--------------<-------------------<------------------<-------------
// 				seq[StartPos-flankRight]........seq[ALtNuc].......seq[endPos+flankLeft]
// 				end <----------------------------------------------------------------start
// 			*/
//
// 			seq = bsatservice.GetRevComplement(bsatservice.GetNucFromGenome(endPos-1, startPos))
// 			refNuc = bsatservice.GetComplement(bsatservice.GetNucFromGenomePos(posToChange))
// 			if flankLeft != 0 {
// 				flankSeqLeft = bsatservice.GetRevComplement(bsatservice.GetNucFromGenome(endPos-1, endPos+flankLeft))
// 			}
// 			if flankRight != 0 {
// 				flankSeqRight = bsatservice.GetNucFromGenome(startPos-1, startPos-flankRight)
// 			}
// 		}
// 		fmt.Println(flankSeqLeft, flankSeqRight, seq, refNuc)
// 		// if bsatstruct.Flag.InfoLocus != "" {
// 		// 	if start != 0 && end != 0 {
// 		// 		switch bsatstruct.Flag.InfoShowAs {
// 		// 		case "", "direct":
// 		// 			seq = bsatservice.GetNucFromGenome(start-1, end)
// 		// 			direction = "in forward direction"
// 		// 		case "gene":
// 		// 			seq = GetGeneSeq(bsatstruct.Flag.InfoLocus)
// 		// 			direction = "in gene direction"
// 		// 			if bsatstruct.Flag.InfoReplace != "" {
// 		// 				posRefAlt := strings.Split(bsatstruct.Flag.InfoReplace, ":")
// 		// 				locPosInGene, _ := strconv.Atoi(posRefAlt[0])
// 		// 				if locPosInGene != 0 && locPosInGene < len(seq) {
// 		// 					splitSeq := strings.Split(seq, "")
// 		// 					if strings.ToUpper(splitSeq[locPosInGene-1]) == strings.ToUpper(posRefAlt[1]) {
// 		// 						if bsatstruct.Flag.InfoIUPAc == false {
// 		//
// 		// 							splitSeq[locPosInGene-1] = fmt.Sprintf("%v", strings.ToUpper(posRefAlt[2]))
// 		// 						} else {
// 		// 							locUIPAC := bsatservice.GetIUPAC(fmt.Sprintf("%v", strings.ToUpper(posRefAlt[2])))
// 		// 							splitSeq[locPosInGene-1] = fmt.Sprintf("\n%v\n", locUIPAC)
// 		// 						}
// 		//
// 		// 					} else {
// 		// 						splitSeq[locPosInGene-1] = fmt.Sprintf("%v", strings.ToUpper(posRefAlt[2]))
// 		// 						fmt.Println("Replaced reference nucleotide does not match to reference!")
// 		// 					}
// 		//
// 		// 					seq = strings.Join(splitSeq, "")
// 		//
// 		// 					if bsatstruct.Flag.InfoFlankRight != 0 || bsatstruct.Flag.InfoFlankLeft != 0 {
// 		// 						seq = fmt.Sprintf("%v%v%v", flankSeqLeft, seq, flankSeqRight)
// 		// 					}
// 		//
// 		// 				}
// 		// 			}
// 		//
// 		// 		case "fasta":
// 		// 			asFasta = true
// 		// 			if bsatstruct.Flag.InfoFlankLeft != 0 {
// 		// 				// flankSeqLeft = bsatservice.GetNucFromGenome((start-bsatstruct.Flag.InfoFlankLeft)-1, end)
// 		// 				start = start - bsatstruct.Flag.InfoFlankLeft
// 		// 				flankInfo.WriteString(fmt.Sprintf("+%v_left ", bsatstruct.Flag.InfoFlankLeft))
// 		// 				// fmt.Println("*", flankSeqLeft, "LEFT")
// 		// 			}
// 		// 			if bsatstruct.Flag.InfoFlankRight != 0 {
// 		// 				flankSeqRight = bsatservice.GetNucFromGenome(end, end+bsatstruct.Flag.InfoFlankRight)
// 		// 				// fmt.Println("*", flankSeqRight, "RIGHT")
// 		// 				flankInfo.WriteString(fmt.Sprintf("+%v_right ", bsatstruct.Flag.InfoFlankRight))
// 		// 			}
// 		//
// 		// 			seq = bsatservice.GetNucFromGenome(start-1, end)
// 		// 			if bsatstruct.Flag.InfoMask != "" {
// 		// 				posRefAlt := strings.Split(bsatstruct.Flag.InfoMask, ":")
// 		// 				locPosInGene, _ := strconv.Atoi(posRefAlt[0])
// 		// 				if locPosInGene != 0 && locPosInGene < len(seq) {
// 		// 					splitSeq := strings.Split(seq, "")
// 		// 					if strings.ToUpper(splitSeq[locPosInGene-1]) == strings.ToUpper(posRefAlt[1]) {
// 		// 						if bsatstruct.Flag.InfoIUPAc == false {
// 		//
// 		// 							splitSeq[locPosInGene-1] = fmt.Sprintf("[%v/%v]", posRefAlt[1], posRefAlt[2])
// 		// 						} else {
// 		// 							locUIPAC := bsatservice.GetIUPAC(fmt.Sprintf("%v%v", posRefAlt[1], posRefAlt[2]))
// 		// 							splitSeq[locPosInGene-1] = fmt.Sprintf("\n%v\n", locUIPAC)
// 		// 						}
// 		//
// 		// 					} else {
// 		// 						splitSeq[locPosInGene-1] = fmt.Sprintf("[%v/%v]", posRefAlt[1], posRefAlt[2])
// 		// 						fmt.Println("Masked reference nucleotide does not match to reference!")
// 		// 					}
// 		//
// 		// 					seq = strings.Join(splitSeq, "")
// 		//
// 		// 					if bsatstruct.Flag.InfoFlankRight != 0 || bsatstruct.Flag.InfoFlankLeft != 0 {
// 		// 						seq = fmt.Sprintf("%v%v%v", flankSeqLeft, seq, flankSeqRight)
// 		// 					}
// 		//
// 		// 				}
// 		// 			}
// 		// 			if bsatstruct.Flag.InfoReplace != "" {
// 		// 				posRefAlt := strings.Split(bsatstruct.Flag.InfoReplace, ":")
// 		// 				locPosInGene, _ := strconv.Atoi(posRefAlt[0])
// 		// 				if locPosInGene != 0 && locPosInGene < len(seq) {
// 		// 					splitSeq := strings.Split(seq, "")
// 		// 					if strings.ToUpper(splitSeq[locPosInGene-1]) == strings.ToUpper(posRefAlt[1]) {
// 		// 						if bsatstruct.Flag.InfoIUPAc == false {
// 		//
// 		// 							splitSeq[locPosInGene-1] = fmt.Sprintf("%v", strings.ToUpper(posRefAlt[2]))
// 		// 						} else {
// 		// 							locUIPAC := bsatservice.GetIUPAC(fmt.Sprintf("%v", strings.ToUpper(posRefAlt[2])))
// 		// 							splitSeq[locPosInGene-1] = fmt.Sprintf("\n%v\n", locUIPAC)
// 		// 						}
// 		//
// 		// 					} else {
// 		// 						splitSeq[locPosInGene-1] = fmt.Sprintf("%v", strings.ToUpper(posRefAlt[2]))
// 		// 						fmt.Println("Replaced reference nucleotide does not match to reference!")
// 		// 					}
// 		//
// 		// 					seq = strings.Join(splitSeq, "")
// 		//
// 		// 					if bsatstruct.Flag.InfoFlankRight != 0 || bsatstruct.Flag.InfoFlankLeft != 0 {
// 		// 						seq = fmt.Sprintf("%v%v%v", flankSeqLeft, seq, flankSeqRight)
// 		// 					}
// 		//
// 		// 				}
// 		// 			}
// 		// 			direction = "in forward direction"
// 		// 		}
// 		//
// 		// 		prod := bsatservice.GetProdByName(bsatstruct.Flag.InfoLocus)
// 		// 		note := bsatservice.GetNoteByName(bsatstruct.Flag.InfoLocus)
// 		// 		goa := bsatservice.GetGOAByName(bsatstruct.Flag.InfoLocus)
// 		//
// 		// 		// вывод последовательности кодонов заданной ключом --codons=Start:End
// 		// 		if bsatstruct.Flag.InfoCodons != "" {
// 		//
// 		// 			seqByCodons := make([]string, (len(seq)/3)+1)
// 		// 			codonNbr := 1
// 		// 			var codonBuffer strings.Builder
// 		// 			for _, c := range seq {
// 		// 				codonBuffer.WriteString(string(c))
// 		// 				if codonBuffer.Len() == 3 {
// 		//
// 		// 					seqByCodons[codonNbr] = codonBuffer.String()
// 		//
// 		// 					codonBuffer.Reset()
// 		// 					codonNbr = codonNbr + 1
// 		// 				}
// 		// 			}
// 		// 			codonsCoords := strings.Split(bsatstruct.Flag.InfoCodons, ":")
// 		// 			startCodon, _ := strconv.Atoi(codonsCoords[0])
// 		// 			endCodon, _ := strconv.Atoi(codonsCoords[1])
// 		// 			for i := startCodon; i <= endCodon; i++ {
// 		// 				if bsatstruct.Flag.GbVerbose == false {
// 		// 					codonBuffer.WriteString(seqByCodons[i])
// 		// 				} else {
// 		//
// 		// 					codonBuffer.WriteString(fmt.Sprintf("[%v]%v", i, seqByCodons[i]))
// 		// 				}
// 		//
// 		// 			}
// 		//
// 		// 			fmt.Printf(">%v %v [%v:%v codons %v]\n%v\n", bsatstruct.Flag.InfoLocus, prod, startCodon, endCodon, direction, codonBuffer.String())
// 		//
// 		// 		} else {
// 		//
// 		// 			if asFasta == false {
// 		// 				gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)
// 		// 				fmt.Printf(
// 		// 					">%v %v (%v-%v %v)\n%v\n-----------------\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n-----------------\n%v\ngoa:%v\n",
// 		// 					bsatstruct.Flag.InfoLocus,
// 		// 					prod,
// 		// 					start, end, direction, seq, gc, gc1, gc2, gc3, len(seq), note, goa)
// 		// 			} else {
// 		//
// 		// 				fmt.Printf(">%v %v (%v-%v %v) %v \n%v\ngoa:%v\n", bsatstruct.Flag.InfoLocus, prod, start, end, direction, flankInfo.String(), seq, goa)
// 		// 			}
// 		//
// 		// 		}
// 		// 	} else {
// 		// 		fmt.Println(colorRed, bsatstruct.Flag.InfoLocus, " not found!")
// 		// 	}
// 		// } else if bsatstruct.Flag.InfoRanges != "" {
// 		// 	coords := strings.Split(bsatstruct.Flag.InfoRanges, ":")
// 		//
// 		// 	start, _ := strconv.Atoi(coords[0])
// 		// 	end, _ := strconv.Atoi(coords[1])
// 		// 	locname, _ := bsatservice.GetGeneNameByPos(start, end)
// 		// 	if start != 0 && end != 0 {
// 		// 		seq := bsatservice.GetNucFromGenome(start-1, end)
// 		// 		// gc, gc1, gc2, gc3 := codon.GcCodonCalc(seq)
// 		// 		// prod := getProductByName(*infoLocus)
// 		// 		// fmt.Printf("%v (%v-%v)\n%v\ngc:%v gc1:%v gc2:%v gc3:%v Len:%v\n", prod, start, end, seq, gc, gc1, gc2, gc3, len(seq))
// 		// 		// fmt.Println(seq, gc, gc1, gc2, gc3, prod)
// 		// 		fmt.Printf(">%v|%v|%v (%v:%v)\n%v\n", bsatstruct.GenomeDescription.Strain, bsatstruct.GenomeDescription.Version, locname, start, end, seq)
// 		//
// 		// 	}
// 		//
// 		// }
// 	}
// }
