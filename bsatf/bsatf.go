package bsatf

import (
	"bsatoolMod/src/amino"
	"bsatoolMod/src/bsatseq"
	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bsatoolMod/src/bsatvcf"
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

/*
Модуль для работы с загрузкой файлов
*/

func LoadExclGenes(file string) map[int]int {
	var (
		exGenes    = make(map[int]int)
		exLocus    = regexp.MustCompile(`(^\w+)`)
		start, end int
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

		scanTxt := scanner.Text()

		for _, exGene := range exLocus.FindAllStringSubmatch(scanTxt, -1) {
			start, end = bsatseq.GetGenePosByName(exGene[0])

			exGenes[start] = end

		}
	}
	return exGenes
}

func LoadExclSNP(file string) map[int]int {
	var (
		exSNPs = make(map[int]int)
		snpPos int
	)

	locSNPcheck := ReadSNPFromFile(file)
	for _, val := range locSNPcheck {
		snpPos, _ = strconv.Atoi(val.APos)
		exSNPs[snpPos] = 1

	}
	return exSNPs
}

func ReadSNPFromFile(f string) []bsatstruct.SnpCheckInfo {
	var (
		parsedSNP    []bsatstruct.SnpCheckInfo
		snpCheckChan = make(chan bsatstruct.SnpCheckInfo)
	)
	file, err := os.Open(f)
	if err != nil {
		log.Fatal(f, err.Error())

	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(f, err.Error())
		}
	}(file)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		if strings.HasPrefix(scanner.Text(), "#") == false || !strings.HasPrefix(scanner.Text(), "//") == false || !strings.HasPrefix(
			scanner.Text(), ";") == false {
			go func() {

				snpCheckChan <- bsatseq.ValidateSNP(scanner.Text())

			}()
			snpCheck := <-snpCheckChan

			parsedSNP = append(parsedSNP, snpCheck)

		}
	}

	return parsedSNP
}

func CheckSNPfromFile(f string, verbose bool, web bool, useRule bool) {
	/*
	   tLMN    = "LMN"
	   tPMLN   = "PMLN"
	   tLSAAN  = "LSAAN"
	   tLLAAN  = "LLAAN"
	*/
	mapOfVCF := make(map[string][]bsatstruct.SnpInfo)
	mapofSNP := make(map[string][]string)

	files := &bsatstruct.ListOfFiles

	locSNPcheck := ReadSNPFromFile(f)

	isNotFound := make(map[string]int)

	snpFound := make(map[string]map[string]int)

	snpFoundCount := make(map[string]int)

	type snpCheckValues struct {
		Locus, NucInPosCoding, RefAAShort, AltAAShort, Tag, Hash string
		APos, CodonNbrInG                                        int
	}

	var (
		chkSNP                   []bsatstruct.CheckSNP
		mutationsLits            []string
		th                       int
		captionMutationList      []string
		checkRes                 = make(map[string][]snpCheckValues)
		mutationRef, mutationAlt string
		mutationPos              int
		snpHashMap               = make(map[uint64]string)
		snpHashList              []uint64
		snpHashPerGenome         = make(map[string][]uint64)
		snpHash                  uint64
		snpHashPerGenomeFound    = make(map[string]map[uint64]int)
		snpHashPerGenomeCount    = make(map[uint64]int)
		snpPosFound              = make(map[string]map[int]int)
		snpPNAFound              = make(map[string]map[int]int)
		// rulesArr                 []bsatstruct.RulesInfo
	)

	if bsatstruct.Flag.StatTH != 0 {
		th = bsatstruct.Flag.StatTH
	} else {
		th = 1
	}

	if useRule == true {
		bsatstruct.RulesArr = CheckRuleFromFile(bsatstruct.Flag.StatInRule)
	}

	for i, file := range *files {
		if verbose == true {
			fmt.Printf("Reading %v (%v:%v)%v\r", file, i+1, len(*files), strings.Repeat(" ", 60))
		}

		snpFound[file] = make(map[string]int)

		snpHashPerGenomeFound[file] = make(map[uint64]int)
		snpPosFound[file] = make(map[int]int)
		snpPNAFound[file] = make(map[int]int)

		qSNP := &bsatvcf.VcfQuery{File: file, OutChan: make(chan bsatvcf.VcfInfoQuery)}
		go qSNP.Request()

		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo

	}

	for fname, snps := range bsatstruct.SnpCacheMap {
		mapOfVCF[fname] = snps
	}

	for file, snp := range mapOfVCF {
		if verbose == true {
			fmt.Printf("Working on %v                     \r", file)
		}

		for _, val := range locSNPcheck {

			CodonNbrInG, _ := strconv.Atoi(val.CodonNbrInG)
			lGpoS, _ := strconv.Atoi(val.PosInGene)
			lAPos, _ := strconv.Atoi(val.APos)

			for _, snpFromFile := range snp {

				switch val.TypeOf {
				// "LMN"   //locus:Mutation:NAME
				// "PMLN"  // position:Mutation:locus:NAME
				// "PMN"   // position:Mutation:NAME 			4326236_C>T	ETO_ethA_Gly413Asp
				// "LSAAN" //locus:shortAA:codon:shortAA:name
				// "LLAAN" // locus:longAA:codon:longAA:name
				// "LCN"   //locus:codon:name  -> Rv0667:codon491:RpoB_codon491
				// "PNA" // position~name~antibotic
				// GLSAAAB // gene(locus)_ShortAACodonShortAA_AB
				case bsatstruct.TGLSAAAB:
					if strings.ToUpper(snpFromFile.Locus) == strings.ToUpper(val.Locus) || strings.ToUpper(snpFromFile.Name) == strings.ToUpper(val.Locus) {

						if snpFromFile.CodonNbrInG == CodonNbrInG && snpFromFile.AltAAShort == val.AASalt {
							mapofSNP[file] = append(
								mapofSNP[file], fmt.Sprintf(
									"%v[%v:%v_%v>%v(%v%v%v)]%v", val.Locus, snpFromFile.Locus,
									snpFromFile.PosInGene,
									strings.ToUpper(snpFromFile.NucInGenomeRef), strings.ToUpper(snpFromFile.NucInGenomeAlt), snpFromFile.RefAAShort,
									snpFromFile.CodonNbrInG, snpFromFile.AltAAShort, val.AB))

							mutationsLits = append(
								mutationsLits, fmt.Sprintf(
									"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.PosInGene,
									strings.ToUpper(snpFromFile.NucInGenomeRef), strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Locus))
							snpFound[file][fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.PosInGene, strings.ToUpper(snpFromFile.NucInGenomeRef),
								strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Locus)] = 1
							snpFoundCount[fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.PosInGene, strings.ToUpper(snpFromFile.NucInGenomeRef),
								strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Locus)] = snpFoundCount[fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus,
								snpFromFile.PosInGene,
								strings.ToUpper(snpFromFile.NucInGenomeRef), strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Locus)] + 1
							snpPNAFound[file][lAPos] = 1
						}
					}

				case bsatstruct.TPNA:
					if lAPos == snpFromFile.APos {
						if snpPNAFound[file][lAPos] != 1 {
							mapofSNP[file] = append(
								mapofSNP[file], fmt.Sprintf(
									"%v[%v:%v_%v>%v]%v", val.Name, snpFromFile.Locus, lAPos,
									strings.ToUpper(snpFromFile.NucInGenomeRef), strings.ToUpper(snpFromFile.NucInGenomeAlt), val.AB))

							mutationsLits = append(
								mutationsLits, fmt.Sprintf(
									"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos,
									strings.ToUpper(snpFromFile.NucInGenomeRef), strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Name))
							snpFound[file][fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(snpFromFile.NucInGenomeRef),
								strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Name)] = 1
							snpFoundCount[fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(snpFromFile.NucInGenomeRef),
								strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Name)] = snpFoundCount[fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos,
								strings.ToUpper(snpFromFile.NucInGenomeRef), strings.ToUpper(snpFromFile.NucInGenomeAlt), val.Name)] + 1
							snpPNAFound[file][lAPos] = 1

						}

					}
				case bsatstruct.TLMN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) && lGpoS == snpFromFile.PosInGene {
						mapofSNP[file] = append(
							mapofSNP[file],
							fmt.Sprintf("%v[%v:%v%v>%v]", val.Name, val.Locus, val.PosInGene, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))

						mutationsLits = append(
							mutationsLits,
							fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1
						snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
							val.Name)] = snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] + 1
					}
				case bsatstruct.TPMLN:
					// "PMLN"  // position:Mutation:locus:NAME
					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) {
						mapofSNP[file] = append(
							mapofSNP[file],
							fmt.Sprintf("%v[%v:%v_%v>%v]", val.Name, snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt)))

						mutationsLits = append(
							mutationsLits,
							fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] = 1
						snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
							val.Name)] = snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] + 1
					}
				case bsatstruct.TLSAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&
						val.AASalt == snpFromFile.AltAAShort {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v%v]", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))
						mutationsLits = append(
							mutationsLits,
							fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
							val.Name)] = 1 // buffer.WriteString(fmt.Sprintf("%v_%v:%v%v%v\t", val.Name, val.Locus, val.AASref, CodonNbrInG, val.AASalt))
						snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
							val.Name)] = snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] + 1

					}
				case bsatstruct.TLLAAN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG &&
						val.AALalt == snpFromFile.AltAA {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:%v%v%v]", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))
						mutationsLits = append(
							mutationsLits,
							fmt.Sprintf("%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name))
						snpFound[file][fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
							val.Name)] = 1 // buffer.WriteString(fmt.Sprintf("%v_%v:%v%v%v\t", val.Name, val.Locus, val.AALref, CodonNbrInG, val.AALalt))
						snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
							val.Name)] = snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, lAPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt), val.Name)] + 1

					}
				case bsatstruct.TLCN:

					if strings.ToUpper(val.Locus) == strings.ToUpper(snpFromFile.Locus) && CodonNbrInG == snpFromFile.CodonNbrInG {
						mapofSNP[file] = append(mapofSNP[file], fmt.Sprintf("%v[%v:codon%v]", val.Name, val.Locus, CodonNbrInG))
						mutationsLits = append(
							mutationsLits, fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.APos,
								strings.ToUpper(snpFromFile.NucInPosCoding),
								strings.ToUpper(snpFromFile.Alt), val.Name))
						snpFound[file][fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPosCoding),
							strings.ToUpper(snpFromFile.Alt),
							val.Name)] = 1 // buffer.WriteString(fmt.Sprintf("%v_%v:codon%v\t", val.Name, val.Locus, CodonNbrInG))
						snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPosCoding),
							strings.ToUpper(snpFromFile.Alt), val.Name)] = snpFoundCount[fmt.Sprintf(
							"%v:%v_%v>%v_%v", snpFromFile.Locus, snpFromFile.APos,
							strings.ToUpper(snpFromFile.NucInPosCoding),
							strings.ToUpper(snpFromFile.Alt), val.Name)] + 1
					}
				case bsatstruct.TPMN:

					// "PMN"   // position:Mutation:NAME
					// 497491_C>T   N:RD105
					// go run bsatool.go stat --db test_core -a check -i drugs5.txt --gene

					if bsatstruct.Flag.StatAsGene == true && snpFromFile.TypeOf == "CDS" {
						mutationRef = snpFromFile.NucInPosCoding
						mutationAlt = snpFromFile.Alt
						mutationPos = snpFromFile.PosInGene
					} else {
						mutationRef = snpFromFile.NucInGenomeRef
						mutationAlt = snpFromFile.NucInGenomeAlt
						mutationPos = lAPos
					}

					snpHash = bsatservice.GetHash(
						fmt.Sprintf(
							"%v%v%v", mutationPos, strings.ToUpper(mutationRef),
							strings.ToUpper(mutationAlt)))

					if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) == strings.ToUpper(snpFromFile.Alt) || lAPos == snpFromFile.APos && strings.
						ToUpper(bsatseq.MakeComplementSeq(val.Alt)) == strings.ToUpper(snpFromFile.Alt) && snpFromFile.Direction == "r" {

						if snpPosFound[file][mutationPos] == 0 {
							mapofSNP[file] = append(
								mapofSNP[file], fmt.Sprintf(
									"%v[%v:%v_%v>%v]", val.Name, snpFromFile.Locus, mutationPos,
									strings.ToUpper(mutationRef),
									strings.ToUpper(mutationAlt)))
							mutationsLits = append(
								mutationsLits, fmt.Sprintf(
									"%v:%v_%v>%v_%v", snpFromFile.Locus, mutationPos, strings.ToUpper(mutationRef),
									strings.ToUpper(mutationAlt), val.Name))
							snpFound[file][fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, mutationPos, strings.ToUpper(snpFromFile.NucInPosCoding),
								strings.ToUpper(val.Alt), val.Name)] = 1
							snpFoundCount[fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, mutationPos, strings.ToUpper(val.Ref), strings.ToUpper(val.Alt),
								val.Name)] = snpFoundCount[fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, mutationPos, strings.ToUpper(val.Ref),
								strings.ToUpper(val.Alt), val.Name)] + 1
							// ---------------------------------------------
							snpHashMap[snpHash] = fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, mutationPos, strings.ToUpper(mutationRef),
								strings.ToUpper(mutationAlt), val.Name)
							snpHashList = append(
								snpHashList, bsatservice.GetHash(
									fmt.Sprintf(
										"%v%v%v", mutationPos, strings.ToUpper(mutationRef),
										strings.ToUpper(mutationAlt))))
							snpHashPerGenome[file] = append(snpHashPerGenome[file], snpHash)
							snpHashPerGenomeFound[file][snpHash] = 1
							snpHashPerGenomeCount[snpHash] = snpHashPerGenomeCount[snpHash] + 1

							// ---------------------------------------------
							snpPosFound[file][mutationPos] = 1
						}

					} else if lAPos == snpFromFile.APos && strings.ToUpper(val.Alt) != strings.ToUpper(snpFromFile.Alt) || lAPos == snpFromFile.
						APos && strings.ToUpper(bsatseq.MakeComplementSeq(val.Alt)) != strings.ToUpper(snpFromFile.Alt) && snpFromFile.Direction == "r" {

						if snpPosFound[file][mutationPos] == 0 {
							mapofSNP[file] = append(
								mapofSNP[file], fmt.Sprintf(
									"%v[%v:%v_%v>%v*]", fmt.Sprintf(
										"%v<-->%v%v%v", val.Name, snpFromFile.RefAA,
										snpFromFile.CodonNbrInG, snpFromFile.AltAA),
									snpFromFile.Locus,
									mutationPos,
									strings.ToUpper(mutationRef),
									strings.ToUpper(mutationAlt)))

							// ---------------------------------------------

							snpHashMap[snpHash] = fmt.Sprintf(
								"%v:%v_%v>%v_%v", snpFromFile.Locus, mutationPos, strings.ToUpper(mutationRef),
								strings.ToUpper(mutationAlt), fmt.Sprintf(
									"%v[%v:%v_%v>%v*]", fmt.Sprintf(
										"%v<-->%v%v%v", val.Name, snpFromFile.RefAA,
										snpFromFile.CodonNbrInG, snpFromFile.AltAA),
									snpFromFile.Locus,
									mutationPos,
									strings.ToUpper(mutationRef),
									strings.ToUpper(mutationAlt)))
							snpHashList = append(snpHashList, snpHash)
							snpHashPerGenome[file] = append(snpHashPerGenome[file], snpHash)

							// ---------------------------------------------

							snpPosFound[file][mutationPos] = 1
						}

					}
				case bsatstruct.TSEN:
					if val.StartRange < snpFromFile.APos && val.EndRange > snpFromFile.APos {
						mapofSNP[file] = append(
							mapofSNP[file], fmt.Sprintf(
								"%v[%v:%v_%v>%v]", snpFromFile.Locus, snpFromFile.Locus, snpFromFile.APos,
								strings.ToUpper(snpFromFile.NucInPosCoding),
								strings.ToUpper(snpFromFile.Alt)))
						mutationsLits = append(
							mutationsLits, fmt.Sprintf(
								"%v:%v%v%v[%v%v%v](%v)", snpFromFile.Locus, strings.ToUpper(snpFromFile.NucInPosCoding),
								snpFromFile.APos, strings.ToUpper(snpFromFile.Alt), snpFromFile.RefAAShort, snpFromFile.CodonNbrInG, snpFromFile.AltAAShort,
								val.Tag))

						// mutationsLits = append(mutationsLits, fmt.Sprintf("%v:%v_%v>%v_%v", fmt.Sprintf("%v[%v:%v_%v>%v]", snpFromFile.Locus, snpFromFile.Locus,
						// 	snpFromFile.APos, strings.ToUpper(snpFromFile.NucInPosCoding),	strings.ToUpper(snpFromFile.Alt))))
						snpFound[file][fmt.Sprintf(
							"%v:%v%v%v[%v%v%v](%v)", snpFromFile.Locus, strings.ToUpper(snpFromFile.NucInPosCoding),
							snpFromFile.APos, strings.ToUpper(snpFromFile.Alt), snpFromFile.RefAAShort, snpFromFile.CodonNbrInG, snpFromFile.AltAAShort,
							val.Tag)] = 1
						snpFoundCount[fmt.Sprintf(
							"%v:%v%v%v[%v%v%v](%v)", snpFromFile.Locus, strings.ToUpper(snpFromFile.NucInPosCoding),
							snpFromFile.APos, strings.ToUpper(snpFromFile.Alt), snpFromFile.RefAAShort, snpFromFile.CodonNbrInG, snpFromFile.AltAAShort,
							val.Tag)] = snpFoundCount[fmt.Sprintf(
							"%v:%v%v%v[%v%v%v](%v)", snpFromFile.Locus, strings.ToUpper(snpFromFile.NucInPosCoding),
							snpFromFile.APos, strings.ToUpper(snpFromFile.Alt), snpFromFile.RefAAShort, snpFromFile.CodonNbrInG, snpFromFile.AltAAShort,
							val.Tag)] + 1
						checkRes[file] = append(
							checkRes[file], snpCheckValues{Locus: snpFromFile.Locus, NucInPosCoding: snpFromFile.NucInPosCoding,
								APos: snpFromFile.APos, RefAAShort: snpFromFile.RefAAShort, AltAAShort: snpFromFile.AltAAShort,
								CodonNbrInG: snpFromFile.CodonNbrInG, Tag: val.Tag})

					}

				}
			}

			isNotFound[file] = len(mapofSNP[file])

		}

		if isNotFound[file] == 0 && len(snpHashPerGenome[file]) == 0 {
			mapofSNP[file] = append(mapofSNP[file], "NO_SNPS_FOUND")

		}

	}

	if len(snpHashList) != 0 {

		for fname, snpUINs := range snpHashPerGenome {

			snpUINs = bsatservice.RmIntDoubles(snpUINs)
			for _, snpUIN := range snpUINs {

				mutationsLits = append(mutationsLits, snpHashMap[snpUIN])
				snpFound[fname][snpHashMap[snpUIN]] = 1
				snpFoundCount[snpHashMap[snpUIN]] = snpFoundCount[snpHashMap[snpUIN]] + 1
			}

		}

	}

	if useRule == true {

		found := make(map[string][]string)
		foundName := make(map[string][]string)

		for fname, val := range mapofSNP {

			for _, tag := range bsatstruct.RulesArr {

				for i := 0; i <= len(tag.Variants)-1; i++ {
					if strings.Contains(strings.ToUpper(strings.Join(val, ",")), tag.Variants[i]) == true {

						found[fname] = bsatservice.AppendIfMiss(found[fname], strings.ToUpper(tag.Variants[i]))

						if len(found[fname]) == tag.Lenght {

							foundName[fname] = bsatservice.AppendIfMiss(foundName[fname], strings.ToUpper(tag.Name))

						}

					}
				}

			}

		}

		for key, val := range mapofSNP {

			sort.Slice(
				found[key], func(i, j int) bool {
					return found[key][i] < found[key][j]
				})

			chkSNP = append(
				chkSNP, bsatstruct.CheckSNP{FileName: key, FoundSNP: strings.Join(val, ","), ExactWithRule: strings.Join(found[key], ","),
					RuleNames: strings.Join(
						foundName[key], ",")})
		}

	} else {
		for key, val := range mapofSNP {
			chkSNP = append(chkSNP, bsatstruct.CheckSNP{FileName: key, FoundSNP: strings.Join(val, ",")})

		}
	}

	if bsatstruct.Flag.StatAsTable == false {

		if web == true {

			bsatservice.PrintSNPFile(chkSNP, bsatstruct.Flag.GbPort)
		} else {

			for _, key := range chkSNP {

				fmt.Printf("%v\t%v\t%v\t%v\n", key.FileName, key.FoundSNP, key.ExactWithRule, key.RuleNames)

			}

		}
	} else {

		// удаляем дубликаты

		mutationsLits = bsatservice.RmStrDoubles(mutationsLits)

		for i := 0; i < len(mutationsLits); i++ {

			if snpFoundCount[mutationsLits[i]] >= th {

				captionMutationList = append(captionMutationList, mutationsLits[i])
			}
		}

		if len(captionMutationList) != 0 {
			// добавляем file в начало captionMutationList
			captionMutationList = append([]string{"file"}, captionMutationList...)
			results := make([]string, len(captionMutationList))
			// выводим заголовки имен мутаций
			fmt.Println(strings.Join(captionMutationList, "\t"))
			for _, file := range *files {
				for i := 0; i < len(captionMutationList); i++ {
					results[0] = file

					if snpFound[file][captionMutationList[i]] == 1 {

						results[i] = "1"
					} else if snpFound[file][captionMutationList[i]] == 0 {
						results[i] = "0"

					}

				}
				fmt.Println(strings.Join(results, "\t"))
			}
		} else {
			fmt.Println("No SNP found in range")
		}

	}

}

func LoadExclRegion(file string) map[int]int {
	var (
		exGenes     = make(map[int]int)
		coordsType1 = regexp.MustCompile(`^\w+\W+(\d+)\W+(\d+)`)

		start, end int
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

		scanTxt := scanner.Text()

		for _, exGene := range coordsType1.FindAllStringSubmatch(scanTxt, -1) {
			start, _ = strconv.Atoi(exGene[1])
			end, _ = strconv.Atoi(exGene[2])

			exGenes[start] = end

		}
	}
	return exGenes
}

func ExclGeneByPos(file string) map[int]int {

	var (
		exGenes = make(map[int]int)
		snpPos  int
	)

	locSNPcheck := ReadSNPFromFile(file)
	for _, val := range locSNPcheck {
		snpPos, _ = strconv.Atoi(val.APos)

		start, end := bsatseq.GetGeneNameCoordsByApos(snpPos)
		if start != 0 && end != 0 {

			exGenes[start] = end
		}

	}
	// fmt.Println(exSNPs)
	return exGenes
}

func GetListVcf() []string {
	/*
		возвращает список VCF файлов в папке в виде массива
	*/

	var (
		fList []string
	)
	files, err := ioutil.ReadDir(".")
	if err != nil {
		log.Fatal(files, err.Error())
	}
	for _, f := range files {

		if strings.ToUpper(path.Ext(f.Name())) == ".VCF" {
			fList = append(fList, f.Name())
		}

	}

	sort.Strings(fList)

	return fList
}

func CheckRuleFromFile(file string) (rules []bsatstruct.RulesInfo) {

	var (
		res   []bsatstruct.RulesInfo
		rinfo bsatstruct.RulesInfo
	)
	rRule := regexp.MustCompile(`^(.*)\b=\b(\w+)$`)

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
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		if !strings.HasPrefix(scanner.Text(), "#") {
			for _, rrule := range rRule.FindAllStringSubmatch(scanner.Text(), -1) {

				convert := strings.Split(strings.ToUpper(rrule[1]), ",")

				sort.Slice(
					convert, func(i, j int) bool {
						return convert[i] < convert[j]
					})
				rinfo = bsatstruct.RulesInfo{Name: rrule[2], Variants: convert, Lenght: len(convert)}

			}
			res = append(res, rinfo)
		}
	}

	return res

}

func AnnotateFromList(file string, genes []bsatstruct.GeneInfo) {

	var (
		posAlt = regexp.MustCompile(`^(\d+)\W+(\w)\W+(\w)`)
		snp    bsatstruct.SnpInfo
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
		for _, findPosAlt := range posAlt.FindAllStringSubmatch(scanner.Text(), -1) {

			if len(findPosAlt) == 4 {

				for z := 0; z < len(genes); z++ {

					apos, _ := strconv.Atoi(findPosAlt[1])

					alt := findPosAlt[3]
					ref := findPosAlt[2]
					if apos >= genes[z].Start && apos <= genes[z].End {
						qSnpInfo := &bsatvcf.SnpInfoQuery{OutChan: make(chan bsatstruct.SnpInfo), Apos: apos, G: genes[z], Alt: alt, Index: false}
						go qSnpInfo.Request()
						snp = <-qSnpInfo.OutChan
						if strings.ToUpper(ref) == strings.ToUpper(snp.NucInGenomeRef) {

							fmt.Printf(
								"%v\t%v\t%v%v>%v\t%v/%v\t%v%v%v\t%v\t%v\n", snp.Locus, snp.APos, snp.PosInGene, snp.NucInPosCoding, snp.Alt, snp.RefCodon,
								snp.AltCodon, snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation, snp.Product)
						} else {

							fmt.Printf(
								"%v\t%v\t%v[%v]>%v\t%v/%v\t%v%v%v\t%v\t%v\n", snp.Locus, snp.APos, snp.PosInGene, snp.NucInPosCoding,
								snp.Alt,
								snp.RefCodon,
								snp.AltCodon, snp.RefAAShort, snp.CodonNbrInG, snp.AltAAShort, snp.Mutation, snp.Product)
						}

					}
				}

			}

		}
	}

}

func CalcComplexIdxFromFile(file string) {

	f, err := os.Open(file) // открываем файл

	if err != nil {
		log.Fatal(file, err.Error())

	}

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		scanTxt := scanner.Text()
		if strings.Contains(scanTxt, "\t") {

			aaArray := strings.Split(scanTxt, "\t")
			if len(aaArray[0]) == 3 && len(aaArray[1]) == 3 {
				idxVal, idxRes := amino.GetComplexIndex(amino.GetShortNameAA(aaArray[0]), amino.GetShortNameAA(aaArray[1]), false)
				fmt.Printf("%v\t%v\t%v\t%v\n", aaArray[0], aaArray[1], idxVal, idxRes)
			} else if len(aaArray[0]) == 1 && len(aaArray[1]) == 1 {
				idxVal, idxRes := amino.GetComplexIndex(aaArray[0], aaArray[1], false)
				fmt.Printf("%v\t%v\t%v\t%v\n", aaArray[0], aaArray[1], idxVal, idxRes)
			}

		}
	}
	defer func(f *os.File) {
		err := f.Close()
		if err != nil {
			log.Fatal(file, err.Error())
		}
	}(f)
}

func FilterVcf(f string, dpFilter int) {
	var (
		vcf = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)

		validateVCF    = regexp.MustCompile(`(##fileformat)=VCF`)
		vcfSpecialInfo = regexp.MustCompile(`(##.*)`)

		vcfValid bool
	)

	file, err := os.Open(f)
	if err != nil {
		log.Fatal(f, err.Error())
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(f, err.Error())
		}
	}(file)

	fOut, err := os.Create("tmp")
	if err != nil {
		log.Fatal("Cannot create file", err.Error())
	}
	defer func(fOut *os.File) {
		err := fOut.Close()
		if err != nil {
			log.Fatal(f, err.Error())
		}
	}(fOut)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {

		for _, vcfvalid := range validateVCF.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfvalid[0] == "##fileformat=VCF" {
				vcfValid = true
			}

		}

		for _, scpecInfo := range vcfSpecialInfo.FindAllStringSubmatch(scanner.Text(), -1) {

			_, _ = fmt.Fprintln(fOut, scpecInfo[0])

		}
		if vcfValid == false {
			fmt.Printf("\n%v is not VCF file!!! Check it!\n", file.Name())
			break
		}
		if bsatstruct.Flag.GbInDel == true {

		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {

				dp, _ := strconv.Atoi(match[4])

				if dp >= dpFilter {
					_, _ = fmt.Fprintln(fOut, scanner.Text())
				}

			}
		}

	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}
	_ = os.Rename(f, fmt.Sprintf("%v.bak", f))
	_ = os.Rename("tmp", f)

}

func CheckPosList(file string) {
	var (
		posInGene   int
		refNuc, seq string
		seqArr      []string
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

		scanTxt := scanner.Text()
		if pos, err := strconv.Atoi(scanTxt); err == nil {
			locus, _ := bsatseq.GetGeneNameByPos(pos, pos)
			prod := bsatseq.GetProdByName(locus)
			start, end := bsatseq.GetGenePosByName(locus)

			if bsatseq.GetDirectByName(locus) == "f" {
				refNuc = bsatseq.GetNucFromGenomePos(pos)
				posInGene = (pos - start) + 1
				seq = bsatseq.GetNucFromGenome(start, end)
				seqArr = strings.Split(seq, "")
				seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", seqArr[posInGene-1], refNuc)
			} else if bsatseq.GetDirectByName(locus) == "r" {
				refNuc = bsatseq.MakeComplementSeq(bsatseq.GetNucFromGenomePos(pos))
				posInGene = (end - pos) + 1
				seq = bsatseq.MakeRevComplement(bsatseq.GetNucFromGenome(start, end))
				seqArr = strings.Split(seq, "")
				seqArr[posInGene-1] = fmt.Sprintf("\n%v(%v)[%v|%v]\n", pos, posInGene, seqArr[posInGene-1], refNuc)
			}

			seq = strings.Join(seqArr, "")
			fmt.Printf("\n>%v %v %v", locus, prod, seq)

		}
	}
}

func CheckPosListVCF(fileUniq string, makeSeq bool) {

	var (
		posCounter                                                                                                        = make(map[int]int)
		foundPos                                                                                                          = make(map[int]bsatstruct.SnpInfo)
		posInGene                                                                                                         int
		refNuc, altNuc, seq, locus, locusDescr, posfromfile, leftFlankSeq, rightFlankSeq, leftFlankDescr, rightFlankDescr string
		seqArr                                                                                                            []string
		posFromFile                                                                                                       = regexp.MustCompile(`(\d+)`)
	)

	files := &bsatstruct.ListOfFiles

	f, err := os.Open(fileUniq) // открываем файл

	if err != nil {
		log.Fatal(fileUniq, err.Error())

	}
	defer func(f *os.File) {
		err := f.Close()
		if err != nil {
			log.Fatal(fileUniq, err.Error())
		}
	}(f)

	scanner := bufio.NewScanner(f) //  новый сканер

	for scanner.Scan() {

		scanTxt := scanner.Text()
		for _, getPos := range posFromFile.FindAllStringSubmatch(scanTxt, -1) {
			posfromfile = getPos[0]
		}

		if pos, err := strconv.Atoi(posfromfile); err == nil {
			for _, file := range *files {

				snpsChan := make(chan []bsatstruct.SnpInfo)

				go func() {
					snpsChan <- bsatvcf.MakeSnps(file)
				}()
				snps := <-snpsChan

				for _, val := range snps {

					if val.APos == pos {
						posCounter[pos] = posCounter[pos] + 1
						if len(*files) == posCounter[pos] {
							foundPos[pos] = val
						}

					}

				}

			}
		}

	}

	var keys []int
	for key := range foundPos {
		keys = append(keys, key)
	}
	sort.Ints(keys)

	if makeSeq == false {

		for _, key := range keys {

			bsatservice.PrintInConsole(foundPos[key], false)

		}
	} else {
		for _, key := range keys {

			if foundPos[key].TypeOf == "CDS" || foundPos[key].TypeOf == "repeat_region" {

				direction := foundPos[key].Direction
				locus = foundPos[key].Locus
				start, end := bsatseq.GetGenePosByName(locus)
				prod := bsatseq.GetProdByName(locus)

				if direction == "f" {
					if bsatstruct.Flag.StatFlankLeft != 0 {
						leftPos := start - bsatstruct.Flag.StatFlankLeft
						leftFlankSeq = bsatseq.GetNucFromGenome(leftPos, end)
						leftFlankDescr = fmt.Sprintf("+left=%v", bsatstruct.Flag.StatFlankLeft)
					}
					if bsatstruct.Flag.StatFlankRight != 0 {
						rightPos := end + bsatstruct.Flag.StatFlankLeft
						rightFlankSeq = bsatseq.GetNucFromGenome(start, rightPos)
						rightFlankDescr = fmt.Sprintf("+right=%v", bsatstruct.Flag.StatFlankRight)
					}
					refNuc = foundPos[key].NucInPosCoding
					posInGene = (key - start) + 1
					altNuc = foundPos[key].Alt
					// posInGeneT := foundPos[key].PosInGene
					seq = bsatseq.GetNucFromGenome(start, end)
					seqArr = strings.Split(seq, "")
					seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", refNuc, altNuc)
					locusDescr = fmt.Sprintf(
						"%v [apos:%v,posInG:%v,iupac:%v, %v %v] ", locus, key, posInGene, bsatservice.GetIUPAC(fmt.Sprintf("%v%v", refNuc, altNuc)),
						leftFlankDescr,
						rightFlankDescr)
				} else if direction == "r" && bsatstruct.Flag.StatReverse == true {
					if bsatstruct.Flag.StatFlankLeft != 0 {
						rightPos := end + bsatstruct.Flag.StatFlankLeft
						rightFlankSeq = bsatseq.MakeRevComplement(bsatseq.GetNucFromGenome(end, rightPos))
						leftFlankDescr = fmt.Sprintf("+c_left=%v", bsatstruct.Flag.StatFlankLeft)
					}
					if bsatstruct.Flag.StatFlankRight != 0 {
						leftPos := start - bsatstruct.Flag.StatFlankLeft
						leftFlankSeq = bsatseq.MakeRevComplement(bsatseq.GetNucFromGenome(leftPos, start))
						rightFlankDescr = fmt.Sprintf("+c_right=%v", bsatstruct.Flag.StatFlankRight)
					}
					refNuc = foundPos[key].NucInPosCoding
					altNuc = foundPos[key].Alt
					posInGene = (end - key) + 1
					seq = bsatseq.MakeRevComplement(bsatseq.GetNucFromGenome(start, end))
					seqArr = strings.Split(seq, "")

					seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", refNuc, altNuc)
					locusDescr = fmt.Sprintf(
						"%v [apos:%v,posInG:%v,iupac:%v, %v %v ] ", locus, key, posInGene, bsatservice.GetIUPAC(fmt.Sprintf("%v%v", refNuc, altNuc)),
						leftFlankDescr,
						rightFlankDescr)

				} else if direction == "r" && bsatstruct.Flag.StatReverse == false {
					if bsatstruct.Flag.StatFlankLeft != 0 {
						leftPos := start - bsatstruct.Flag.StatFlankLeft
						leftFlankSeq = bsatseq.GetNucFromGenome(leftPos, end)
						leftFlankDescr = fmt.Sprintf("+left=%v", bsatstruct.Flag.StatFlankLeft)
					}
					if bsatstruct.Flag.StatFlankRight != 0 {
						rightPos := end + bsatstruct.Flag.StatFlankLeft
						rightFlankSeq = bsatseq.GetNucFromGenome(start, rightPos)
						rightFlankDescr = fmt.Sprintf("+right=%v", bsatstruct.Flag.StatFlankRight)
					}
					refNuc = bsatseq.MakeRevComplement(foundPos[key].NucInPosCoding)
					posInGene = (end - key) + 1
					altNuc = bsatseq.MakeRevComplement(foundPos[key].Alt)
					// posInGeneT := foundPos[key].PosInGene
					seq = bsatseq.GetNucFromGenome(start, end)
					seqArr = strings.Split(seq, "")
					seqArr[posInGene-1] = fmt.Sprintf("\n[%v/%v]\n", refNuc, altNuc)
					locusDescr = fmt.Sprintf(
						"%v [apos:%v,posInG:%v,iupac:%v, %v %v] ", locus, key, posInGene, bsatservice.GetIUPAC(fmt.Sprintf("%v%v", refNuc, altNuc)),
						leftFlankDescr,
						rightFlankDescr)
				}

				seq = strings.Join(seqArr, "")
				if bsatstruct.Flag.StatFlankLeft != 0 {
					seq = fmt.Sprintf("%v%v", leftFlankSeq, seq)
					leftFlankSeq = ""
				}
				if bsatstruct.Flag.StatFlankRight != 0 {
					seq = fmt.Sprintf("%v%v", seq, rightFlankSeq)
					rightFlankSeq = ""
				}
				fmt.Printf("\n>%v %v \n%v\n", locusDescr, prod, seq)
			}
		}
	}
}

func MakeSeqFromSNPListFile(f string) {

	var (
		snpCheckChan = make(chan bsatstruct.SnpCheckInfo)
		seq          string
	)
	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			log.Fatal(f, err.Error())
		}
	}(file)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {

		if !strings.HasPrefix(scanner.Text(), "#") {

			go func() {

				snpCheckChan <- bsatseq.ValidateSNP(scanner.Text())
			}()
			snpCheck := <-snpCheckChan
			if snpCheck.Locus == "" {
				aposInt, _ := strconv.Atoi(snpCheck.APos)
				locName, _ := bsatseq.GetGeneNameByPos(aposInt, aposInt)
				snpCheck.Locus = locName
			}

			if snpCheck.Locus != "" {

				dir := bsatseq.GetDirectByName(snpCheck.Locus)

				seq = bsatseq.GetGeneSeq(snpCheck.Locus)

				prod := bsatseq.GetProdByName(snpCheck.Locus)
				seqArr := strings.Split(seq, "")
				posInGene, _ := strconv.Atoi(snpCheck.PosInGene)
				for i := 0; i < len(seqArr); i++ {
					if bsatstruct.Flag.GbVerbose == true {

						seqArr[i] = fmt.Sprintf("%v->%v,", i+1, seqArr[i])
					}
					if i+1 == posInGene {
						if dir == "f" {
							seqArr[i] = fmt.Sprintf("[%v>%v]", strings.ToUpper(seqArr[i]), snpCheck.Alt)
						} else if dir == "r" {
							seqArr[i] = fmt.Sprintf("[%v>%v]", strings.ToUpper(seqArr[i]), bsatseq.MakeComplementSeq(snpCheck.Alt))
						}
					}
				}

				altSeq := strings.Join(seqArr, "")

				fmt.Printf(">%v_%v_%v>%v (%v) %v\n%v\n", snpCheck.Locus, snpCheck.PosInGene, snpCheck.Ref, snpCheck.Alt, prod, snpCheck.Name, altSeq)

			}
		}
	}

}

func LoadFileNames(file string) (FileList []string) {

	f, err := os.Open(file) // открываем файл

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

		vcfFile := scanner.Text()
		// fmt.Println(vcfFile)
		if _, err := os.Stat(vcfFile); os.IsNotExist(err) {
			fmt.Printf("%v is not exist\n", vcfFile)

		} else {
			FileList = append(FileList, vcfFile)
		}
	}

	return FileList
}
