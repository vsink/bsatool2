package bsatseq

import (
	"bsatoolMod/src/bsatservice"
	"bsatoolMod/src/bsatstruct"
	"bufio"
	"log"
	"os"
	"regexp"
	"sort"
	"strconv"
	// "bsatoolMod/src/bsatvcf"
	// "bsatoolMod/src/bsatservice"
	codon "bsatoolMod/src/codonPkg"
	"fmt"
	// "log"
	// "math/rand"
	// "os"
	// "sort"
	// "strconv"
	"strings"
	// "time"
)

func MakeComplementSeq(sequence string) string {

	seqArr := strings.Split(sequence, "")

	for i := 0; i < len(seqArr); i++ {
		switch seqArr[i] {
		case "a":
			seqArr[i] = "t"
		case "c":
			seqArr[i] = "g"
		case "t":
			seqArr[i] = "a"
		case "g":
			seqArr[i] = "c"
		case "A":
			seqArr[i] = "T"
		case "C":
			seqArr[i] = "G"
		case "T":
			seqArr[i] = "A"
		case "G":
			seqArr[i] = "C"

		}
	}
	complStr := strings.Join(seqArr, "")

	return complStr

}

func GetDirectByName(locus string) string {
	var direction string

	for _, g := range bsatstruct.FullGenesInfo {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			direction = g.Direction

			break
		}
	}

	return direction
}

func GetGeneNameByPos(start, end int) (string, int) {
	var (
		locName string
		locLen  int
	)
	locLen = (end - start) + 1

	for _, g := range bsatstruct.FullGenesInfo {
		if start >= g.Start && end <= g.End {
			locName = g.Locus

		}

	}
	if locName == "" {
		locName = "IGR"
	}
	return locName, locLen
}

func GetProdByPos(start, end int) (string, string) {
	var prod, note string

	for _, g := range bsatstruct.FullGenesInfo {

		if start >= g.Start && end <= g.End {

			prod = g.Product
			note = g.Note

		}

	}

	return prod, note

}

func GetGeneNameCoordsByApos(pos int) (int, int) {
	var (
		gStart, gEnd int
	)

	for _, g := range bsatstruct.FullGenesInfo {
		if pos >= g.Start && pos <= g.End {
			gStart = g.Start
			gEnd = g.End

		}

	}

	return gStart, gEnd
}

func GetGenePosByName(locus string) (int, int) {
	var start, end int

	g := bsatstruct.GenePositions[locus]

	if g.Start != 0 {
		start = g.Start
		end = g.End
	}

	return start, end
}

func GetGeneSeq(locus string) string {
	var seq string
	start, end := GetGenePosByName(locus)
	if start != 0 && end != 0 {
		dir := GetDirectByName(locus)
		if dir == "f" {
			seq = GetNucFromGenome(start-1, end)
		} else if dir == "r" {
			seq = MakeRevComplement(GetNucFromGenome(start-1, end))

		}
	} else {
		fmt.Println("Locus not found [getGeneSequence func]", locus)
	}
	return seq
}

func GetNucFromGenome(start int, end int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var (
		result string
		slice  []string
	)
	slice = bsatstruct.GenomeSequence[start:end]
	result = strings.Join(slice, "")
	return result
}

func MakeRevComplement(sequence string) string {

	seqArr := strings.Split(sequence, "")

	for i := 0; i < len(seqArr); i++ {
		switch seqArr[i] {
		case "a":
			seqArr[i] = "t"
		case "c":
			seqArr[i] = "g"
		case "t":
			seqArr[i] = "a"
		case "g":
			seqArr[i] = "c"
		case "A":
			seqArr[i] = "T"
		case "C":
			seqArr[i] = "G"
		case "T":
			seqArr[i] = "A"
		case "G":
			seqArr[i] = "C"

		}
	}
	complStr := strings.Join(seqArr, "")

	runes := []rune(complStr)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}

	return string(runes)

}

func GetProdByName(locus string) string {
	var prod string

	for _, g := range bsatstruct.FullGenesInfo {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			prod = g.Product
			break
		}
	}

	return prod
}

func MakeAltString(locus string, positions []bsatstruct.AllPositionsInGene) string {
	// var lStart, lEnd int
	var (
		seqSplit []string
		seq      string
	)

	seq = GetGeneSeq(locus)

	for _, nuc := range seq {
		seqSplit = append(seqSplit, string(nuc))

	}

	for _, val := range positions {
		// fmt.Println(val.pos, val.ref, ">", val.alt, seqSplit[val.pos-1], ">", val.alt)
		if len(seqSplit) != 0 {
			seqSplit[val.Pos-1] = val.Alt
		}

	}

	return strings.Join(seqSplit, "")

}

func GetNucFromGenomePos(pos int) string {
	/*
		возвращает последовательность нуклеотидов из генома:

	*/
	var (
		result string
		slice  []string
	)
	slice = bsatstruct.GenomeSequence[pos-1 : pos]
	result = strings.Join(slice, "")
	return result

}

func MakeAltStringByPos(start int, end int, positions []bsatstruct.AllPositionsInGene) string {
	// var lStart, lEnd int
	var (
		seqSplit []string
		seq      string
	)

	seq = GetNucFromGenome(start, end)

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

func GetCoordRange(start, end int) {

	var coordArray []int

	coordArray = append(coordArray, start)

	for i := start; i <= end; i++ {

		for _, g := range bsatstruct.FullGenesInfo {

			if g.Start == i {

				coordArray = append(coordArray, g.Start)

			} else if g.End == i {

				coordArray = append(coordArray, g.End)

			}

		}
	}
	coordArray = append(coordArray, end)
	sort.Ints(coordArray)

	var res []bsatstruct.RangePosInfo

	for _, val := range coordArray {

		res = append(res, GetSequenceRange(fmt.Sprintf("%v:%v", val, val), bsatstruct.Flag.GbNoSeq))

	}

	var last string
	for _, val := range res {
		if val.GeneName != last {
			fmt.Println(val.GeneName, val.Prod)
		}
		last = val.GeneName

	}

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

func GetSequenceRange(PosRange string, noseq bool) bsatstruct.RangePosInfo {

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
			gname, _ := GetGeneNameByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			prod, _ := GetProdByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			result = bsatstruct.RangePosInfo{Start: bsatstruct.ResStart + 1, End: bsatstruct.ResEnd, GeneName: gname, Prod: prod}

		} else if noseq == false {
			seq := GetNucFromGenome(bsatstruct.ResStart-1, bsatstruct.ResEnd)
			gname, _ := GetGeneNameByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			prod, _ := GetProdByPos(bsatstruct.ResStart, bsatstruct.ResEnd)
			result = bsatstruct.RangePosInfo{Start: bsatstruct.ResStart + 1, End: bsatstruct.ResEnd, GeneName: gname, Len: len(seq), Seq: seq, Prod: prod}

		}

	} else {
		fmt.Println("Range positions is not valid")
	}

	return result

}

func GetRangeFromFile(file string, verbose bool, noseq bool, genome string) []bsatstruct.RangePosInfo {
	var (
		rangeParser = regexp.MustCompile(`\b(\d+)\b\W+\b(\d+)\b`)

		posRange, seq string
		gc            float64
		unsorted      []bsatstruct.RangeArray
		result        []bsatstruct.RangePosInfo

		j, k int
	)

	if genome == "" {
		genome = bsatstruct.GenomeDescription.Strain
	}
	f, err := os.Open(file)
	if err != nil {
		log.Fatal(err)
	}
	defer func(f *os.File) {
		err := f.Close()
		if err != nil {
			log.Fatal(err)
		}
	}(f)

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {

		posRange = scanner.Text()

		posRange = strings.Replace(posRange, ".", "", -1)

		for _, val := range rangeParser.FindAllStringSubmatch(posRange, -1) {

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
			if bsatstruct.ResStart != 0 && bsatstruct.ResEnd != 0 {

				currRange := bsatstruct.RangeArray{Start: bsatstruct.ResStart, End: bsatstruct.ResEnd}
				unsorted = append(unsorted, currRange)
				j++
				if verbose == true {
					fmt.Printf("Processed %v ranges...\r", j)
				}

			}

		}

	}

	sort.Slice(
		unsorted, func(i, j int) bool {
			return unsorted[i].Start < unsorted[j].Start
		})

	encountered := map[int]bool{}
	doubles := map[int]int{}
	var found []bsatstruct.RangeArray

	for v := range unsorted {
		if encountered[unsorted[v].Start+unsorted[v].End] == true {
			// Do not add duplicate.
			doubles[unsorted[v].Start+unsorted[v].End] = doubles[unsorted[v].Start+unsorted[v].End] + 1
		} else {
			// Record this element as an encountered element.
			encountered[unsorted[v].Start+unsorted[v].End] = true
			// Append to result slice.
			found = append(found, unsorted[v])
		}
	}

	for _, val := range found {
		if noseq == true {
			seq = ""
		} else {
			seq = GetNucFromGenome(val.Start-1, val.End)
			if len(seq) > 3 {
				gc, _, _, _ = codon.GcCodonCalc(seq)
			}
		}
		gname, _ := GetGeneNameByPos(val.Start, val.End)
		prod, note := GetProdByPos(val.Start, val.End)
		fixedProd := strings.Replace(prod, " + ", " ", -1)
		gcRes, _ := strconv.ParseFloat(fmt.Sprintf("%.2f", gc), 64)

		res := bsatstruct.RangePosInfo{Start: val.Start, End: val.End, GeneName: gname, Prod: fixedProd, Len: val.End - val.Start + 1, Seq: seq,
			Doubles: doubles[val.Start+val.End] + 1, Note: note, GC: gcRes, Genome: genome}

		result = append(result, res)
		k++
		if verbose == true {
			fmt.Printf("Annotated %v ranges from %v %v\r", k, len(found)-1, strings.Repeat(" ", 10))
		}
	}

	return result
}

func ChkPosExists(s []bsatstruct.AllPositionsInGene, pos int, alt string) bool {
	for _, a := range s {
		if a.Pos == pos && a.Alt == alt {
			return true
		}
	}
	return false
}

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

func GetGenomeMap() []bsatstruct.GenomeMapInfo {

	var (
		gmap    []bsatstruct.GenomeMapInfo
		gmapval bsatstruct.GenomeMapInfo
	)
	for _, g := range bsatstruct.FullGenesInfo {
		gmapval = bsatstruct.GenomeMapInfo{Start: g.Start, End: g.End, Locus: g.Locus, TypeOf: g.TypeOf, Product: g.Product}
		gmap = append(gmap, gmapval)
	}

	return gmap
}

func GetGOAByName(locus string) string {
	var goa string

	for _, g := range bsatstruct.FullGenesInfo {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			goa = g.GOA
			break
		}
	}

	return goa
}

func GetInterGen(pos int) {
	var (
		igensS []int
		igensE []int
	)
	for i, g := range bsatstruct.FullGenesInfo {

		igensS = append(igensS, g.Start)
		igensE = append(igensE, g.End)

		if i >= 1 && i < len(bsatstruct.FullGenesInfo)-1 {

			fmt.Printf("igen:s%v e:%v\n", igensS[i], igensE[i-1])
		}
		fmt.Printf("l:%v s:%v e:%v d:%v %v\n", g.Locus, g.Start, g.End, g.Direction, i)

	}
}

func GetNoteByName(locus string) string {
	var note string

	for _, g := range bsatstruct.FullGenesInfo {
		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) {
			note = g.Note
			break
		}
	}

	return note
}

func Nuc2IntCode(nuc string) string {

	var (
		locNuc = strings.ToUpper(nuc)
		res    string
	)

	switch locNuc {

	}

	switch locNuc {

	case "A":
		res = "0"

	case "T":
		res = "1"
	case "G":
		res = "2"
	case "C":

		res = "3"
	}
	return res
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
					// if gstrain == "" {
					// 	strain = val2[1]
					// } else {
					// 	strain = gstrain
					// }
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
						// fmt.Println("yes", strain, end, first, pos)
						positions = append(positions, pos)
						sequence = append(sequence, nuc)

					} else {

						if len(positions) != 0 {
							// fmt.Println(positions[0], positions[len(positions)-1], strain)
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
				// if gstrain == "" {
				// 	strain = val[1]
				// } else {
				// 	strain = gstrain
				// }
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

				// if count == pos {

				// }

				// fmt.Println(strain, first, end, fmt.Sprintf("%v(%v)L%v",gname, dp, end-first))

				// fmt.Printf("s:%v p:%v s:%v\n", first, pos, pos+1)
				//
				// if (pos == first+1) == false {
				// 	fmt.Printf("s:%v p:%v s:%v\n", first, pos, pos+1)
				// }

			}
			// first = pos

		}

		// if err := scanner.Err(); err != nil {
		// 	log.Fatal(err)
		// }

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

		// toCircos(allGenesVal)
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

// func reverse(data []string) []string {
// 	for i := 0; i < len(data)/2; i++ {
// 		j := len(data) - i - 1
// 		data[i], data[j] = data[j], data[i]
// 	}
// 	return data
// }

// func reverseString(s string) string {
// 	rns := []rune(s) // convert to rune
// 	for i, j := 0, len(rns)-1; i < j; i, j = i+1, j-1 {
//
// 		// swap the letters of the string,
// 		// like first with last and so on.
// 		rns[i], rns[j] = rns[j], rns[i]
// 	}
//
// 	// return the reversed string.
// 	return string(rns)
// }

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

func GetInfoByLocus() {
	var (
		seq, direction              string
		asFasta                     bool
		flankSeqLeft, flankSeqRight string
		flankInfo                   strings.Builder
		colorRed                    = "\x1b[31;1m"
	)
	start, end := GetGenePosByName(bsatstruct.Flag.InfoLocus)
	if bsatstruct.Flag.InfoLocus != "" {
		if start != 0 && end != 0 {
			switch bsatstruct.Flag.InfoShowAs {
			case "", "direct":
				seq = GetNucFromGenome(start-1, end)
				direction = "in forward direction"
			case "gene":
				seq = GetGeneSeq(bsatstruct.Flag.InfoLocus)
				direction = "in gene direction"
			case "fasta":
				asFasta = true
				if bsatstruct.Flag.InfoFlankLeft != 0 {
					flankSeqLeft = GetNucFromGenome((start-bsatstruct.Flag.InfoFlankLeft)-1, end)
					flankInfo.WriteString(fmt.Sprintf("+%v_left ", bsatstruct.Flag.InfoFlankLeft))
				}
				if bsatstruct.Flag.InfoFlankRight != 0 {
					flankSeqRight = GetNucFromGenome(end, end+bsatstruct.Flag.InfoFlankRight)
					flankInfo.WriteString(fmt.Sprintf("+%v_right ", bsatstruct.Flag.InfoFlankRight))
				}

				seq = GetNucFromGenome(start-1, end)
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

			prod := GetProdByName(bsatstruct.Flag.InfoLocus)
			note := GetNoteByName(bsatstruct.Flag.InfoLocus)
			goa := GetGOAByName(bsatstruct.Flag.InfoLocus)

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
		locname, _ := GetGeneNameByPos(start, end)
		if start != 0 && end != 0 {
			seq := GetNucFromGenome(start-1, end)
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
			prod = GetProdByName(g.Locus)
			gc, gc1, gc2, gc3 = codon.GcCodonCalc(seq)
			fmt.Printf("%v:%v\t%v\t%v\t%.2f\t%.2f\t%.2f\t%.2f\n", start, end, g.Locus, prod, gc, gc1, gc2, gc3)

		}
	}
}
