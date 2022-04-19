package bsatvcf

import (
	"bsatoolMod/src/amino"
	"bsatoolMod/src/bsatseq"
	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bufio"
	"fmt"
	"github.com/pkg/browser"
	"html/template"
	"log"
	"math/rand"
	"net/http"
	"os"
	"regexp"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"
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

func (q *VcfQuery) Request() {
	q.OutChan <- VcfInfoQuery{File: q.File, SnpInfo: ParserVCF(q.File, false, bsatstruct.Flag.GbDP, bsatstruct.FullGenesInfo)}

}
func (q *VcfQuery) GetOutChan() <-chan VcfInfoQuery {
	return q.OutChan
}

func (q *SnpInfoQuery) Request() {
	// q.OutChan = make(chan snpInfo)
	q.OutChan <- GetSnpInfo(q.Apos, q.G, q.Alt, q.Index)

}

func (q *SnpInfoQuery) GetOutChan() <-chan bsatstruct.SnpInfo {
	return q.OutChan
}

func GetSnpInfo(apos int, g bsatstruct.GeneInfo, alt string, flgTang bool) bsatstruct.SnpInfo {

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
		codonVal = bsatseq.GetNucFromGenome((posInGene+lStart)-2, ((posInGene+lStart)-1)+2)
	} else if posInCodonG == 1 {
		codonVal = bsatseq.GetNucFromGenome((posInGene+lStart)-3, ((posInGene+lStart)-1)+1)
	} else if posInCodonG == 2 {
		codonVal = bsatseq.GetNucFromGenome((posInGene+lStart)-4, (posInGene+lStart)-1)
	}

	/*
		ревертируем кодон для комплементарных генов

	*/
	nucG := bsatseq.GetNucFromGenomePos((posInGene + lStart) - 1)
	genomePos := nucG
	genomeAlt := alt
	typeOf = g.TypeOf

	if g.Direction == "r" {

		codonVal = bsatseq.MakeRevComplement(codonVal)
		posInGene = CPosInGene
		codonNbrInG = CCodonNbrInG
		alt = bsatseq.MakeComplementSeq(alt)
		nucG = bsatseq.MakeComplementSeq(nucG)
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
							go qSnpInfo.Request()
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
						go qSnpInfo.Request()
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

func MakeSnps(fname string) (snps []bsatstruct.SnpInfo) {
	qSNP := &VcfQuery{File: fname, OutChan: make(chan VcfInfoQuery)}
	go qSNP.Request()

	snpRes := <-qSNP.OutChan
	snps = snpRes.SnpInfo
	return snps
}

func MakeSeq(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) []bsatstruct.SeqInfo {

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
		altPercent     int
		altMinMax      []string
		altMin, altMax int
	)

	if bsatstruct.Flag.AnnSeqLen != 0 {
		nbrOfSNP = bsatstruct.Flag.AnnSeqLen
	}

	files := &bsatstruct.ListOfFiles

	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := &VcfQuery{File: file, OutChan: make(chan VcfInfoQuery), Print: verbose}
		go qSNP.Request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
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

						// fmt.Println(val.APos)

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
			if len(dpMAP[key]) >= bsatstruct.Flag.AnnMinPosNbr {
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
			if altPercent < altMin {
				altPercent = altMin
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

	// fmt.Println(SelectedPos)
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
				refBuffer.WriteString(bsatseq.GetNucFromGenomePos(allpos))
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

					buffer.WriteString(bsatseq.GetNucFromGenomePos(allpos))

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

func NucWebServer(port string, exGenes map[int]int, exSNPs map[int]int) {
	/*

	 */

	seq := MakeSeq(bsatstruct.NCFlg, bsatstruct.Flag.GbVerbose, bsatstruct.Flag.AnnMakeSeqRef, exGenes, exSNPs, bsatstruct.Flag.GbRandomize)
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

func GetAltPos(start int, end int, vcfFile string) []bsatstruct.AllPositionsInGene {
	var (
		snps         []bsatstruct.SnpInfo
		altPositions []bsatstruct.AllPositionsInGene
	)

	snpsChan := make(chan []bsatstruct.SnpInfo)

	go func() {
		snpsChan <- MakeSnps(vcfFile)
	}()
	snps = <-snpsChan

	for _, val := range snps {

		if start >= val.Start && end <= val.End && val.TypeOf == "CDS" {
			altPositions = append(altPositions, bsatstruct.AllPositionsInGene{Pos: val.PosInGene, Alt: val.Alt, Ref: val.NucInPosCoding, Locus: val.Locus})

		}
	}

	return altPositions
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
// 		snpsChan <- MakeSnps(vcfFile)
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

func Pos2Seq(file string, vcfFile string) bsatstruct.AltStringResult {
	var (
		coords = regexp.MustCompile(`^(\S+)\W+(\d+)\W+(\d+)`)

		start, end          int
		locus, prod, altseq string
		altPostitions       []bsatstruct.AllPositionsInGene
		result              bsatstruct.AltStringResult
	)

	f, err := os.Open(file) // открываем файл

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

			locus, _ = bsatseq.GetGeneNameByPos(start, end)
			prod, _ = bsatseq.GetProdByPos(start, end)
			altPostitions = GetAltPos(start, end, vcfFile)

			altseq = bsatseq.MakeAltStringByPos(start, end, altPostitions)

			result = bsatstruct.AltStringResult{Start: start, End: end, Locus: locus, AltSeq: altseq, Prod: prod, VcfFile: vcfFile}

		}

	}
	return result
}

func Snp2SeqByLoc(locus string, vcfFile string) bsatstruct.AltStringResult {
	var (
		start, end    int
		prod, altseq  string
		altPostitions []bsatstruct.AllPositionsInGene
		result        bsatstruct.AltStringResult
	)

	start, end = bsatseq.GetGenePosByName(locus)
	prod, _ = bsatseq.GetProdByPos(start, end)
	altPostitions = GetAltPos(start, end, vcfFile)

	altseq = bsatseq.MakeAltString(locus, altPostitions)

	result = bsatstruct.AltStringResult{Start: start, End: end, Locus: locus, AltSeq: altseq, Prod: prod, VcfFile: vcfFile}

	return result
}

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
		snpsChan <- MakeSnps(vcfFile)
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
			prod := bsatseq.GetProdByName(allloc)
			fmt.Println(allloc, prod, len(altPositions[allloc]))
		}

	}

}

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
		go qSNP.Request()

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

func GetUniqSNP(snps []bsatstruct.SnpInfo, exGenes map[int]int, exSNPs map[int]int) map[int]int {
	var uniqSNP = make(map[int]int)
	for _, val := range snps {
		if len(exGenes) != 0 {

			for key, value := range exGenes {
				if val.APos >= key && val.APos <= value {

					uniqSNP[val.APos] = 2
					continue
				} else if exSNPs[val.APos] == 1 {

					uniqSNP[val.APos] = 2
					continue

				} else {
					if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
						uniqSNP[val.APos] = 1

					}
				}
			}
		} else {
			if uniqSNP[val.APos] != 2 && uniqSNP[val.APos] != 1 {
				uniqSNP[val.APos] = 1

			}

		}

	}
	return uniqSNP
}

func MakeAlign(verbose bool, ref bool) {

	type posPerFile struct {
		Fname string
		Pos   []map[int]string
	}

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
		qSNP := VcfQuery{File: file, OutChan: make(chan VcfInfoQuery), Print: verbose}
		go qSNP.Request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo

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
				lineSeq = append(lineSeq, bsatseq.GetNucFromGenomePos(i))
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

func MakeSeqBinary(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool, dp int) {

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
		qSNP := VcfQuery{File: file, OutChan: make(chan VcfInfoQuery), Print: verbose}
		go qSNP.Request()
		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
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

			if len(dpMAP[key]) >= bsatstruct.Flag.AnnMinPosNbr {
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

func MakeSeqNex(typeof string, verbose bool, ref bool, exGenes map[int]int, exSNPs map[int]int, randomize bool) {

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

		// excludedLocus = make(map[int]string)
		// excludedProd  = make(map[int]string)
		// excludedDP    = make(map[int][]int)

		// posCount             = make(map[int]int)
	)

	// files := GetListVcf()

	// queryChan := make(chan vcfInfoQuery)

	if bsatstruct.Flag.AnnSeqLen != 0 {
		nbrOfSNP = bsatstruct.Flag.AnnSeqLen
	}

	files := &bsatstruct.ListOfFiles
	for i, file := range *files {
		aaAltCoords[file] = make(map[int]string)

		// создаем запрос в виде типа vcfQuery, который передается через канал на выход <-qSNP.OutChan
		qSNP := VcfQuery{File: file, OutChan: make(chan VcfInfoQuery), Print: verbose}
		go qSNP.Request()
		// snpCacheFromChan = append(snpCacheFromChan, <-qSNP.OutChan)
		snpRes := <-qSNP.OutChan
		bsatstruct.SnpCacheMap[snpRes.File] = snpRes.SnpInfo
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

						// fmt.Println(val.APos)

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

	// fmt.Println(uniqueSNP)

	for key, value := range uniqueSNP {

		if value == 1 {
			if len(dpMAP[key]) >= bsatstruct.Flag.AnnMinPosNbr {
				AllPos = append(AllPos, key)

			}

		}

	}
	// go process("Working...          ")
	// AllPos = unique(AllPosUnsort)

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
			// rnd := rand.Intn(len(AllPos)-i) + i
			// fmt.Println(AllPos[rnd])

			TempPos = append(TempPos, AllPos[i])
		}
	}

	sort.Ints(TempPos)
	// fmt.Println(posFreq)
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
			// if *annMinPosNbr!=0 &&  {
			//
			// }
			altPercent = (count1 * 100) / len(posFreq[pos])
			altMinMax = strings.Split(bsatstruct.Flag.AnnAltRange, ":")
			if len(altMinMax) != 0 {
				altMin, _ = strconv.Atoi(altMinMax[0])
				altMax, _ = strconv.Atoi(altMinMax[1])
			}
			// fmt.Println(altPercent,altMin,altMax)
			if altPercent >= altMin && altPercent <= altMax {
				SelectedPos = append(SelectedPos, pos)
				if bsatstruct.Flag.GbDebug {
					fmt.Printf("pos:%v isRef:%v,isAlt:%v %v legnth_array: %v AltPerc: %v \n", pos, count0, count1, posFreq[pos], len(posFreq[pos]), altPercent)
					// fmt.Println(*annTest)
				}
			}

		}

	}

	// -pos_file ФЛАГ

	sort.Ints(SelectedPos)
	for fname, snps := range bsatstruct.SnpCacheMap {

		pos := make(map[int]string)
		var buffer strings.Builder
		// var aaSNPpos []int

		buffer.WriteString(fmt.Sprintf(">%v\n", fname))
		switch typeof {
		case bsatstruct.NCFlg:
			for _, val := range snps {
				pos[val.APos] = val.Alt

			}
			nTax = len(bsatstruct.SnpCacheMap)
			nChar = len(SelectedPos)
			for _, allpos := range SelectedPos {
				// posCount[allpos] = posCount[allpos] + 1
				if pos[allpos] != "" {

					// buffer.WriteString(pos[allpos])
					nexusTaxa[fname] = append(nexusTaxa[fname], pos[allpos])

				} else {

					// buffer.WriteString(getNucFromGenomePos(allpos))
					nexusTaxa[fname] = append(nexusTaxa[fname], bsatseq.GetNucFromGenomePos(allpos))
				}

			}

		}

	}
	if len(nexusTaxa) != 0 {
		// var taxNbr int
		// fmt.Println(nexusTaxa)
		// fmt.Printf("#nexus\nBEGIN Taxa;\nDIMENSIONS\nntax=%v;\nTAXLABELS\n", nTax)
		fmt.Printf("#NEXUS\nBegin data;\nDimensions ntax=%v nchar=%v\nFormat datatype=dna missing=N gap=-;\n", nTax, nChar)
		// for key := range nexusTaxa {
		// 	taxNbr++
		// 	fmt.Printf("[%v]\t'%v'\n", taxNbr, key)
		// }
		fmt.Printf("matrix\n")
		for key, val := range nexusTaxa {

			fmt.Printf("%v  \n%v\n", key, strings.Join(val, ""))
		}
		fmt.Printf(";\nEnd;\n")
	}
	// fmt.Println(buffer.Len())
	// fmt.Println(len(buffer.String()))

	// 	case aaFlag:
	//
	// 		for _, allpos := range SelectedPos {
	// 			if aaAltCoords[fname][allpos] == "" {
	// 				buffer.WriteString(strings.ToLower(aaRefCoords[allpos]))
	// 				aaSNPpos = append(aaSNPpos, allpos)
	// 			} else if aaAltCoords[fname][allpos] != "" {
	// 				aaSNPpos = append(aaSNPpos, allpos)
	// 				buffer.WriteString(aaAltCoords[fname][allpos])
	// 			}
	//
	// 		ResSeq = append(ResSeq, seqInfo{Name: fname, Seq: buffer.String(), UsedPositions: aaSNPpos})
	//
	// 	}
	// }

	// return ResSeq
}

func ParserBulkVCF(withFilenames bool) {

	files := &bsatstruct.ListOfFiles
	for i, file := range *files {

		snpsChan := make(chan []bsatstruct.SnpInfo)

		go func() {
			snpsChan <- MakeSnps(file)
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
