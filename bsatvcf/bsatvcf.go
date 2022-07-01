package bsatvcf

import (
	"bsatoolMod/src/amino"
	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
)

type (
	VcfInfoQuery struct {
		File    string
		SnpInfo []bsatstruct.SnpInfo
	}

	VcfQuery struct {
		OutChan chan VcfInfoQuery
		File    string
		Print   bool
	}
	SnpInfoQuery struct {
		OutChan chan bsatstruct.SnpInfo
		Apos    int
		G       bsatstruct.GeneInfo
		Alt     string
		Index   bool
	}
)

func (q *VcfQuery) request() {
	q.OutChan <- VcfInfoQuery{File: q.File, SnpInfo: ParserVCF(q.File, false, bsatstruct.Flag.GbDP, bsatstruct.FullGenesInfo)}

}
func (q *VcfQuery) getOutChan() <-chan VcfInfoQuery {
	return q.OutChan
}

func (q *SnpInfoQuery) request() {
	// q.OutChan = make(chan snpInfo)
	q.OutChan <- getSnpInfo(q.Apos, q.G, q.Alt, q.Index)

}

func (q *SnpInfoQuery) getOutChan() <-chan bsatstruct.SnpInfo {
	return q.OutChan
}

// func GetAltPositionsVCFIGR(start int, end int, vcfFile string) []bsatstruct.AllPositionsInGene {
// 	var (
// 		// allLocuses []string
// 		// locus        string
// 		snps         []bsatstruct.SnpInfo
// 		altPositions []bsatstruct.AllPositionsInGene
// 		// count        = 1
// 	)
//
// 	snpsChan := make(chan []bsatstruct.SnpInfo)
//
// 	go func() {
// 		snpsChan <- GetSNPList(vcfFile)
// 	}()
// 	snps = <-snpsChan
// 	// fmt.Println(locus)
// 	for _, val := range snps {
// 		// fmt.Println(val.Locus, val.PosInGene, val.Alt)
// 		if start >= val.Start && end <= val.End {
// 			altPositions = append(altPositions, bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})
// 			// locus = val.Locus
// 		}
// 	}
//
// 	return altPositions
// }

// func GetHashSNP(snp bsatstruct.SnpInfo) uint64 {
// 	// hasher := sha1.New()
// 	// hasher.Write([]byte(text))
// 	// return hex.EncodeToString(hasher.Sum(nil))
// 	hash, err := hashstructure.Hash(snp, nil)
// 	if err != nil {
// 		panic(err)
// 	}
//
// 	// fmt.Printf("%d", hash)
// 	return hash
// }

func GetShareSNP(verbose bool, igr bool, web bool, files []string) {

	var (
		pos   = map[uint64]int{}
		alt   = map[uint64]bsatstruct.SnpInfo{}
		share []uint64

		snpToWeb     []bsatstruct.SnpInfo
		snpToConsole []bsatstruct.SnpInfo
	)
	upperLimit := len(files)

	for i, file := range files {

		if verbose == true {
			fmt.Printf("Pass: %v from %v \r", i+1, len(files))
		}

		qSNP := VcfQuery{File: file, OutChan: make(chan VcfInfoQuery), Print: verbose}
		go qSNP.request()

		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
	}

	for _, snps := range bsatstruct.SnpCacheMap {

		for _, val := range snps {

			hash := bsatservice.GetHash(fmt.Sprintf("%v%v%v%v%v%v%v", val.Locus, val.APos, val.PosInGene, val.Alt, val.Direction, val.TypeOf, val.Product))

			pos[hash] = pos[hash] + 1

			alt[hash] = val

		}

		for i, lPos := range pos {

			if lPos == upperLimit {

				share = append(share, i)

			}

		}

	}

	for _, sharePos := range share {

		if web == true {
			snpToWeb = append(snpToWeb, alt[sharePos])

		} else {

			snpToConsole = append(snpToConsole, alt[sharePos])
		}
	}

	//  --------------------------------------------
	if web == true && len(snpToWeb) != 0 {
		sort.Slice(
			snpToWeb, func(i, j int) bool {
				return snpToWeb[i].Start < snpToWeb[j].Start
			})
		bsatservice.PrintInWEB(snpToWeb, bsatstruct.Flag.GbPort)
	} else if web == false && len(snpToConsole) != 0 {

		sort.Slice(
			snpToConsole, func(i, j int) bool {
				return snpToConsole[i].Start < snpToConsole[j].Start
			})
		for _, res := range snpToConsole {

			bsatservice.PrintInConsole(res, igr)
		}
	}
	//  ---------------------------------------------

}

// func GetUniqSNP(snps []bsatstruct.SnpInfo, exGenes map[int]int, exSNPs map[int]int) map[int]int {
// 	var uniqSNP = make(map[int]int)
// 	for _, val := range snps {
// 		if len(exGenes) != 0 {
//
// 			for key, value := range exGenes {
// 				if val.APos >= key && val.APos <= value {
//
// 					uniqSNP[val.APos] = 2
// 					continue
// 				} else if exSNPs[val.APos] == 1 {
//
// 					uniqSNP[val.APos] = 2
// 					continue
//
// 				} else {
// 					if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
// 						uniqSNP[val.APos] = 1
//
// 					}
// 				}
// 			}
// 		} else {
// 			if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
// 				uniqSNP[val.APos] = 1
//
// 			}
//
// 		}
//
// 	}
// 	return uniqSNP
// }

func CountSNP(vcfFile string) {

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
		snpsChan <- GetSNPList(vcfFile)
	}()
	snps = <-snpsChan

	for _, val := range snps {
		if val.TypeOf == "CDS" && bsatservice.ChkPosExists(altPositions[val.Locus], val.PosInGene, val.Alt) == false {

			altPositions[val.Locus] = append(
				altPositions[val.Locus], bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})
		}

	}

	for _, allloc := range allLocuses {

		if len(altPositions[allloc]) > 2 {
			prod := bsatservice.GetProdByName(allloc)
			fmt.Println(allloc, prod, len(altPositions[allloc]))
		}

	}

}

func ParserVCF(f string, print bool, dpFilter int, genes []bsatstruct.GeneInfo) []bsatstruct.SnpInfo {
	var (
		vcf         = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*DP=(\d+)`)
		indel       = regexp.MustCompile(`^\S+\s+(\d+)\W+(\w+)\s+(\w+).*INDEL.*DP=(\d+)`)
		validateVCF = regexp.MustCompile(`(##fileformat)=VCF`)
		vcfValid    bool

		//
		snpFromVCF []bsatstruct.SnpInfo
	)
	cpuprofile := bsatstruct.Flag.AnnBench
	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Fatal(err)

		}
		_ = pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	file, err := os.Open(f)
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
		for _, vcfvalid := range validateVCF.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfvalid[0] == "##fileformat=VCF" {
				vcfValid = true
			}

		}
		if vcfValid == false {
			fmt.Printf("\n%v is not VCF file!!! Check it!\n", file.Name())
			break
		}
		if bsatstruct.Flag.GbInDel == true {
			for _, matchindel := range indel.FindAllStringSubmatch(scanner.Text(), -1) {
				if vcfValid == true {
					apos, _ := strconv.Atoi(matchindel[1])
					ref := matchindel[2]
					alt := matchindel[3]
					dp, _ := strconv.Atoi(matchindel[4])

					for z := 0; z < len(genes); z++ {

						if apos >= genes[z].Start && apos <= genes[z].End && dp >= dpFilter {
							if bsatstruct.Flag.GbDebug == true {
								fmt.Printf("apos:%v\tref:%v\talt:%v\tdp: %v\t f:%v\n", apos, ref, alt, dp, f)
							}

							qSnpInfo := SnpInfoQuery{OutChan: make(chan bsatstruct.SnpInfo), Apos: apos, G: genes[z], Alt: alt,
								Index: bsatstruct.Flag.GbIndex}
							go qSnpInfo.request()
							snp := <-qSnpInfo.OutChan

							snp.InterPro = genes[z].InterPro
							snp.PDB = genes[z].PDB
							snp.ProSite = genes[z].ProSite
							snp.DP = dp
							snp.Indel = 1
							snp.Alt = fmt.Sprintf("%v/%v", ref, alt)
							snp.IndelAlt = alt
							snp.IndelRef = ref
							if len(ref) > len(alt) {
								snp.IndelType = fmt.Sprintf("del%v", alt)

							} else if len(ref) < len(alt) {
								snp.IndelType = fmt.Sprintf("ins%v", alt)
							}

							if bsatstruct.Flag.GbDebug == true {
								if strings.ToUpper(snp.NucInPosCoding) == strings.ToUpper(snp.Alt) {
									snp.Product = fmt.Sprintf("%v|!ref=alt[%v]dp=%v", snp.Product, f, dp)
								}

							}

							if strings.ToUpper(snp.NucInPosCoding) != strings.ToUpper(snp.Alt) {
								snpFromVCF = append(snpFromVCF, snp)
							}

							if print == true {

								bsatservice.PrintInConsole(snp, bsatstruct.Flag.GbIGR)

							}

						}

					}
				}
			}
		}

		for _, match := range vcf.FindAllStringSubmatch(scanner.Text(), -1) {

			if vcfValid == true {
				apos, _ := strconv.Atoi(match[1])
				ref := match[2]
				alt := match[3]
				dp, _ := strconv.Atoi(match[4])

				for z := 0; z < len(genes); z++ {

					if apos >= genes[z].Start && apos <= genes[z].End && dp >= dpFilter {

						qSnpInfo := &SnpInfoQuery{OutChan: make(chan bsatstruct.SnpInfo), Apos: apos, G: genes[z], Alt: alt, Index: bsatstruct.Flag.GbIndex}
						go qSnpInfo.request()
						snp := <-qSnpInfo.OutChan

						snp.InterPro = genes[z].InterPro
						snp.PDB = genes[z].PDB
						snp.ProSite = genes[z].ProSite
						snp.DP = dp
						snp.Indel = 0
						if bsatstruct.Flag.GbDebug == true {
							if strings.ToUpper(snp.NucInPosCoding) == strings.ToUpper(snp.Alt) {
								snp.Product = fmt.Sprintf("%v|!ref=alt[%v]dp=%v", snp.Product, f, dp)
							}

						}

						if len(ref) == 1 && len(alt) == 1 && strings.ToUpper(snp.NucInPosCoding) != strings.ToUpper(snp.Alt) {
							snpFromVCF = append(snpFromVCF, snp)
							if print == true {

								bsatservice.PrintInConsole(snp, bsatstruct.Flag.GbIGR)

							}
						}

					}

				}
			}
		}

	}

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	return snpFromVCF

}

func getSnpInfo(apos int, g bsatstruct.GeneInfo, alt string, flgTang bool) bsatstruct.SnpInfo {

	var (
		snp                         bsatstruct.SnpInfo // структура SNP
		codonPositions              []string           // срез для разбивки кодона побуквенно
		altCodonPositions           []string           // срез для разбивки кодона побуквенно альтернативным нуклеотидом
		locReportType, typeOf, titv string
		geneLen                     int
		codonVal                    string // переменная для кодона
		altCodon                    string // аналогично для альтернативного кодона
		mut                         string
		tangIdx                     string
		tangIdxVal                  int
	)
	// var trouble int
	lStart := g.Start // переменная начала гена
	lEnd := g.End
	posInGene := (apos - lStart) + 1             // позиция снипа в гене
	codonNbrInG := ((posInGene - 1) / 3) + 1     // номер кодона=номеру аминокислоты в трансляции
	posInCodonG := (codonNbrInG * 3) - posInGene // позиция в буквы в кодоне (0-первая, 1-средняя, 2-последняя)
	CPosInGene := (lEnd - apos) + 1              // комплементарный ген. позиция в гене
	CCodonNbrInG := ((CPosInGene - 1) / 3) + 1   // комплементарный ген. номер кодона = номеру аминокислоты
	// финт, который делал в snpMiner2, сдвиг на 1 букву. взял оттуда

	if posInCodonG == 2 {
		posInCodonG = 0
	} else if posInCodonG == 0 {
		posInCodonG = 2
	}

	/*
		определяем границы кодона, в зависимости от положения нуклеотида
		0- три буквы справа
		1-одна справа, одна слева
		2-три буквы слева
	*/
	if lStart > lEnd {
		geneLen = (lStart - lEnd) + 1

	} else if lEnd > lStart {
		geneLen = (lEnd - lStart) + 1

	}
	if posInCodonG == 0 {
		codonVal = bsatservice.GetNucFromGenome((posInGene+lStart)-2, ((posInGene+lStart)-1)+2)
	} else if posInCodonG == 1 {
		codonVal = bsatservice.GetNucFromGenome((posInGene+lStart)-3, ((posInGene+lStart)-1)+1)
	} else if posInCodonG == 2 {
		codonVal = bsatservice.GetNucFromGenome((posInGene+lStart)-4, (posInGene+lStart)-1)
	}

	/*
		ревертируем кодон для комплементарных генов

	*/
	nucG := bsatservice.GetNucFromGenomePos((posInGene + lStart) - 1)
	genomePos := nucG
	genomeAlt := alt
	typeOf = g.TypeOf

	if g.Direction == "r" {

		codonVal = bsatservice.GetRevComplement(codonVal)
		posInGene = CPosInGene
		codonNbrInG = CCodonNbrInG
		alt = bsatservice.GetComplement(alt)
		nucG = bsatservice.GetComplement(nucG)
		if posInCodonG == 2 {
			posInCodonG = 0
		} else if posInCodonG == 0 {
			posInCodonG = 2
		}

	}
	//

	codonPositions = strings.Split(codonVal, "")

	codonPositions[posInCodonG] = strings.ToUpper(codonPositions[posInCodonG])

	codonVal = strings.Join(codonPositions, "")

	altCodonPositions = codonPositions

	altCodonPositions[posInCodonG] = alt

	altCodonPositions[posInCodonG] = strings.ToUpper(altCodonPositions[posInCodonG])

	altCodon = strings.Join(altCodonPositions, "")

	aaRef, aaRefShort := amino.Codon2AA(codonVal)
	aaAlt, aaAltShort := amino.Codon2AA(altCodon)
	if aaRefShort == aaAltShort {
		mut = "synonymous"
		tangIdx = "0000"
	} else if aaRefShort != aaAltShort && aaAltShort != "X" {
		mut = "missense"

		// tangIdx = strconv.FormatFloat(amino.GetTangInx(aaRefShort, aaAltShort), 'f', 2, 64)
		tangIdx, tangIdxVal = amino.GetComplexIndex(aaRefShort, aaAltShort, bsatstruct.Flag.GbVerbose)
	} else if aaRefShort != aaAltShort && aaAltShort == "X" {
		mut = "nonsense"
		tangIdx = "0000"
	}

	if flgTang == true {
		locReportType = "T1"
	} else if flgTang == false {
		locReportType = "T0"
	}

	titv = bsatservice.CalcTiTv(nucG, alt)

	snp = bsatstruct.SnpInfo{APos: apos, PosInGene: posInGene, PosInCodonG: posInCodonG,
		RefCodon: codonVal, RefAA: aaRef, NucInPosCoding: strings.ToUpper(nucG), Locus: g.Locus,
		Direction: g.Direction, Name: g.Name, Product: g.Product,
		Start: g.Start, End: g.End, CodonNbrInG: codonNbrInG, AltCodon: altCodon,
		AltAA: aaAlt, RefAAShort: aaRefShort, AltAAShort: aaAltShort,
		Mutation: mut, Tang: tangIdx, TangIdxVal: tangIdxVal, Alt: alt,
		Note: g.Note, ReportType: locReportType, ProteinID: g.ProteinID,
		GeneID: g.GeneID, GOA: g.GOA, GeneLen: geneLen, TiTv: titv, TypeOf: typeOf, NucInGenomeRef: genomePos, NucInGenomeAlt: genomeAlt}

	return snp
}

func GetSNPList(fname string) (snps []bsatstruct.SnpInfo) {
	qSNP := &VcfQuery{File: fname, OutChan: make(chan VcfInfoQuery)}
	go qSNP.request()

	snpRes := <-qSNP.OutChan
	snps = snpRes.SnpInfo
	return snps
}

func GetSNPInfoList(apos int, g bsatstruct.GeneInfo, alt string) (snp bsatstruct.SnpInfo) {
	qSnpInfo := &SnpInfoQuery{OutChan: make(chan bsatstruct.SnpInfo), Apos: apos, G: g, Alt: alt, Index: false}
	go qSnpInfo.request()
	snp = <-qSnpInfo.OutChan
	return snp
}
