package bsatservice

import (
	"bsatoolMod/src/bsatstruct"
	"bsatoolMod/src/codonPkg"
	"bufio"
	"github.com/360EntSecGroup-Skylar/excelize"
	"github.com/pterm/pterm"
	"log"
	"os/exec"
	"regexp"
	"runtime"
	"strconv"

	"fmt"
	"github.com/mitchellh/hashstructure"
	"github.com/pkg/browser"
	"html/template"
	// "log"
	// "math/rand"
	"net/http"
	"os"
	"sort"
	"strings"
	// "time"
)

func CalcTiTv(ref, alt string) string {
	var typeOf string

	if ref != "" && alt != "" {
		lref, lalt := strings.ToUpper(ref), strings.ToUpper(alt)
		if lref == "A" && lalt == "G" || lref == "G" && lalt == "A" || lref == "C" && lalt == "T" || lref == "T" && lalt == "C" {
			typeOf = "Ti"
		} else {
			typeOf = "Tv"
		}
	}
	return typeOf
}

func PrintInConsole(snps bsatstruct.SnpInfo, igr bool) {

	const fullAnnotations = " {{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\") (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T0\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPosCoding}}>{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPosCoding}}>{{.Alt}}\t-\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{end}}"

	const fullAnnotationsWithName = " {{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.FName}}\t{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPosCoding}}>{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (ne .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.APos}}{{.NucInPosCoding}}>{{.Alt}}\t-\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{end}}"

	const cdsAnnotations = "{{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 0))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\")  (eq .Indel 1))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.IndelType}}\t{{.Alt}}\t-\t-\t-\t-\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\")  (eq .Indel 1))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\t{{.IndelType}}\n" +
		"{{end}}"

	const cdsAnnotationsWithName = "{{if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T0\"))}}" +
		"{{.FName}}\t{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Product}}\n" +
		"{{else if (and (eq .TypeOf \"CDS\") (eq .ReportType \"T1\"))}}" +
		"{{.Locus}}\t{{.APos}}\t{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}\t{{.RefCodon}}/{{.AltCodon}}\t{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}\t{{.Mutation}}\t{{.Tang}}({{.TangIdxVal}})\t{{.Product}}\n" +
		"{{end}}"

	t := template.New("report")
	switch bsatstruct.Flag.AnnShowFileName {
	case false:
		if igr == true {
			t, _ = t.Parse(fullAnnotations)
			_ = t.Execute(os.Stdout, snps)

		} else {
			t, _ = t.Parse(cdsAnnotations)
			_ = t.Execute(os.Stdout, snps)
		}
	case true:
		if igr == true {
			t, _ = t.Parse(fullAnnotationsWithName)
			_ = t.Execute(os.Stdout, snps)

		} else {

			t, _ = t.Parse(cdsAnnotationsWithName)
			_ = t.Execute(os.Stdout, snps)
		}
	}

}

func GetHash(str string) uint64 {

	hash, err := hashstructure.Hash(str, nil)
	if err != nil {
		panic(err)
	}

	return hash
}

func RmIntDoubles(elements []uint64) []uint64 {
	var (
		result      []uint64
		encountered = map[uint64]bool{}
	)

	// Create a map of all unique elements.
	for v := range elements {
		encountered[elements[v]] = true
	}

	// Place all keys from the map into a slice.

	for key := range encountered {
		result = append(result, key)
	}
	sort.Slice(result, func(i, j int) bool { return result[i] < result[j] })
	return result
}

func AppendIfMiss(slice []string, val string) []string {
	sort.Slice(
		slice, func(i, j int) bool {
			return slice[i] < slice[j]
		})

	for _, ele := range slice {
		if ele == val {
			return slice
		}
	}

	return append(slice, val)
}

func PrintSNPFile(stat []bsatstruct.CheckSNP, port string) {

	var htmlTemplate = `   <!DOCTYPE html>
			<html>
			<head>
			<meta charset="utf-8">			
			</head>		
			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<body>
			<tr>		
		 
			{{range $element := .}}	
				{{ $length := len .ExactWithRule}}{{ if eq $length 0 }}
			   				<td>{{.FileName}}</td><td>{{.FoundSNP}}</td>		
					{{else}} 
							<td>{{.FileName}}</td><td>{{.FoundSNP}}</td><td>{{.ExactWithRule}}</td>
				{{end}}
				
			</tr>
			{{end}}	
			</body>
			</table>
		
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>			
			`

	t := template.New("t")

	t, err := t.Parse(htmlTemplate)

	if err != nil {
		panic(err)
	}

	_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

	http.HandleFunc(
		"/", func(w http.ResponseWriter, r *http.Request) {
			// err = t.Execute(w, &gInfo)

			err = t.Execute(w, stat)
			if err != nil {
				panic(err)
			}
			go func() {
				defer os.Exit(0)
			}()
		})

	locPort := fmt.Sprintf(":%v", port)
	_ = http.ListenAndServe(locPort, nil)
}

func RmStrDoubles(elements []string) []string {
	var (
		result      []string
		encountered = map[string]bool{}
	)

	// Create a map of all unique elements.
	for v := range elements {
		encountered[elements[v]] = true
	}

	// Place all keys from the map into a slice.

	for key := range encountered {
		result = append(result, key)
	}
	sort.Strings(result)
	return result
}

func PrintInWEB(snps []bsatstruct.SnpInfo, port string) {

	var htmlTitle = `   <!DOCTYPE html>
				<html>

			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<tbody>
			<tr>
			
			{{if eq .NucCore 1}}
			
				<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
				<td>AA</td><td>Type</td><td>Product</td><td><p title="Gene Ontology Annotation (GOA) Database">GOA</p></td><td><p title="The Universal Protein Resource (UniProt) is a comprehensive resource for protein sequence and annotation data">UniProt</p></td><td><p title="InterPro provides functional analysis of proteins by classifying them into families and predicting domains and important sites">InterPro</p></td><td><p title="Protein Data Bank (PDB)">PDB</p></td><td>ProSite</td>
			{{else if eq .NucCore 0}}
				<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
				<td>AA</td><td>Type</td><td>Product</td><td>GeneID</td><td>UniProt</td>
			{{else if  eq .NucCore 2}}
				<td>Locus</td><td>Gene</td><td>Pos.</td><td>Mutation</td><td>Codons</td>
				<td>AA</td><td>Type</td><td>Product</td><td>-</td><td>UniProt</td>
			{{end}}
			
			</tr>
			
`
	var htmlTemplate = `
			
			{{/* 1 */}}
				{{range $element := .}}
				{{/* 2 */}}
			
			{{if .GeneID}}
				<tr>
					
				<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}</td>
				<td>{{.RefCodon}}/{{.AltCodon}}</td><td><p title="{{.RefAA}}{{.CodonNbrInG}}{{.AltAA}}">{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</p></td>
				<td><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}"target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ncbi.nlm.nih.gov/gene/{{.GeneID}}={{.GeneID}}" target="_blank">{{.GeneID}}</a>
				</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score">{{.ProteinID}}</td>
				</tr>
					{{/* 2 */}}
			{{else}}
				<tr>
				
					{{/* 3 */}}
			{{if eq .TypeOf "CDS"}}			
				<td><p title="{{.Note}}">{{.Locus}}</p></td><td>{{.Name}}</td><td>{{.APos}}</td><td>{{.PosInGene}}{{.NucInPosCoding}}>{{.Alt}}</td>
				<td>{{.RefCodon}}/{{.AltCodon}}</td><td><p title="{{.RefAA}}{{.CodonNbrInG}}{{.AltAA}}">{{.RefAAShort}}{{.CodonNbrInG}}{{.AltAAShort}}</p></td>			
					{{/* 4 */}}
			{{if eq .Mutation "missense"}}
					{{if eq .TangIdxVal 4 }}
				<td bgcolor="#CECEF6"><u><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></u></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{else}}
	<td bgcolor="#CECEF6"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{end}}
					{{/* 4 */}}
			{{else if eq .Mutation "nonsense"}}
				<td bgcolor="#F78181"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{/* 4 */}}
			{{else}}
				<td bgcolor="#ddffcc"><p title="Complex Index: {{.Tang}}({{.TangIdxVal}})">{{.Mutation}}</p></td><td><a href="https://www.ncbi.nlm.nih.gov/protein/{{.ProteinID}}" target="_blank"><p title="{{.Note}}">{{.Product}}</p></a></td><td><a href="https://www.ebi.ac.uk/QuickGO/GProtein?ac={{.GOA}}" target="_blank">{{.GOA}}</a>
					{{/* 4 */}}
			{{end}}
				</td><td><a href="http://www.uniprot.org/uniprot/?query={{.ProteinID}}&sort=score" target="_blank">{{.ProteinID}}</td>
				<td>
				{{ range $value := .InterPro }}
   				<a href="http://www.ebi.ac.uk/interpro/entry/{{$value}}" target="_blank">{{$value}}</a>
				{{ end }}
				</td>
				<td>
				{{ range $value := .PDB }}
   				<a href="https://www.rcsb.org/structure/{{$value}}" target="_blank">{{$value}}</a>
				{{ end }}			
				</td>
				<td>
				{{ range $value := .ProSite }}
   				<a href="https://prosite.expasy.org/{{$value}}" target="_blank">{{$value}}</a>
				{{ end }}			
				</td>
				</tr>
				{{/* 3 */}}
			{{end}}
				{{/* 2 */}}
			{{end}}
				{{/* 1 */}}		
			{{end}}
			
			</tbody>
			</table>
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>
			`

	t := template.New("t")
	tt := template.New("t1")

	t, err := t.Parse(htmlTemplate)
	t1, err := tt.Parse(htmlTitle)

	if err != nil {
		panic(err)
	}

	locURL := fmt.Sprintf("localhost:%v/", port)
	_ = browser.OpenURL(locURL)

	http.HandleFunc(
		"/", func(w http.ResponseWriter, r *http.Request) {
			// err = t.Execute(w, &gInfo)

			err = t1.Execute(w, &bsatstruct.GenomeDescription)
			if err != nil {
				panic(err)
			}
			err = t.Execute(w, snps)
			if err != nil {
				panic(err)
			}
			go func() {
				defer os.Exit(0)
			}()
		})

	// locPort := fmt.Sprintf(":%v", port)
	_ = http.ListenAndServe(":8080", nil)
}

// func AnnotatePos(file string) {
// 	cols := regexp.MustCompile(`(\w+)\W+(\d+)\W+(\d+)`)
// 	f, err := os.Open(file)
//
// 	if err != nil {
// 		fmt.Println(err)
//
// 	}
// 	defer f.Close()
// 	fmt.Println("name\tstart\tend\tlocus\tlen\tgc\tseq\tproduct")
// 	scanner := bufio.NewScanner(f)
//
// 	for scanner.Scan() {
// 		for _, colsMatch := range cols.FindAllStringSubmatch(scanner.Text(), -1) {
// 			// fmt.Printf("c1:%v c2:%v\n", colsMatch[1], colsMatch[2])
// 			lStart, _ := strconv.Atoi(colsMatch[2])
// 			lEnd, _ := strconv.Atoi(colsMatch[3])
// 			gname, glen := bsatseq.GetGeneNameByPos(lStart, lEnd)
// 			// prod := getProductByName(gname)
//
// 			prod, _ := bsatseq.GetProdByPos(lStart, lEnd)
// 			seq := bsatseq.GetNucFromGenome(lStart-1, lEnd)
// 			gc, _, _, _ := codon.GcCodonCalc(seq)
// 			fmt.Printf("%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", colsMatch[1], colsMatch[2], colsMatch[3], gname, glen, gc, seq, prod)
//
// 		}
// 	}
//
// }

// func ChkLocus(locus string) int {
// 	var locType int
//
// 	for _, g := range bsatstruct.FullGenesInfo {
// 		if strings.ToUpper(locus) == strings.ToUpper(g.Locus) && g.TypeOf == "CDS" {
// 			locType = 1
// 			break
// 		} else if strings.ToUpper(locus) == strings.ToUpper(g.Locus) && g.TypeOf != "CDS" {
// 			locType = 2
// 			break
//
// 		} else if strings.ToUpper(locus) == strings.ToUpper(g.Locus) && g.TypeOf == "" {
// 			locType = 0
// 			break
//
// 		}
// 	}
//
// 	return locType
// }

func GetIUPAC(nucleotides string) string {

	var iupac string
	locNuc := strings.ToTitle(nucleotides)

	if locNuc == "AC" || locNuc == "CA" {
		iupac = "M"
	} else if locNuc == "AG" || locNuc == "GA" {
		iupac = "R"
	} else if locNuc == "CT" || locNuc == "TC" {
		iupac = "Y"
	} else if locNuc == "GC" || locNuc == "CG" {
		iupac = "S"
	} else if locNuc == "AT" || locNuc == "TA" {
		iupac = "W"
	} else if locNuc == "GT" || locNuc == "TG" {
		iupac = "K"
	}
	return iupac
}

func CompareSlices(slice1 []string, slice2 []string) []string {
	var diffStr []string
	m := map[string]int{}

	for _, s1Val := range slice1 {
		m[s1Val] = 1
	}
	for _, s2Val := range slice2 {
		m[s2Val] = m[s2Val] + 1
	}

	for mKey, mVal := range m {
		if mVal == 1 {
			diffStr = append(diffStr, mKey)
		}
	}

	return diffStr
}

func PrintStatsInWeb(stat []bsatstruct.StatInfo, port string) {
	var htmlTemplate = `   <!DOCTYPE html>
				<html>
				<head>
				<style>
				.col {
				word-wrap: break-word; /* Перенос слов */
				}
				</style>
			</head>
			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<tbody>
			<tr>			{{range $element := .}}
			{{if eq .Perc 100}}				
			<td>{{.Pos}}</td><td>{{.Count}}</td><td bgcolor="#ffe6e6">{{.Perc}}%</td><td>ALL FILES</td><td>NONE</td>
			{{else if eq .Perc 5}}				
			<td>{{.Pos}}</td><td>{{.Count}}</td><td bgcolor="#ddffcc">{{.Perc}}%</td><td>{{.FilesWith}}</td><td>{{.FilesWithout}}</td>
			{{else}}
			<td>{{.Pos}}</td><td>{{.Count}}</td><td>{{.Perc}}%</td><td>{{.FilesWith}}</td><td>{{.FilesWithout}}</td>
			{{end}}
			<td>Pos</td><td>Count</td><td>Percent</td><td>+</td><td>-</td>
			</tr>
			<tr>
			{{range $element := .}}
			{{if eq .Perc 100}}				
			<td>{{.Pos}}</td><td>{{.Count}}</td><td bgcolor="#ffe6e6">{{.Perc}}%</td><td>ALL FILES</td><td>NONE</td>
			{{else if eq .Perc 5}}				
			<td>{{.Pos}}</td><td>{{.Count}}</td><td bgcolor="#ddffcc">{{.Perc}}%</td><td>{{.FilesWith}}</td><td>{{.FilesWithout}}</td>
			{{else}}
			<td>{{.Pos}}</td><td>{{.Count}}</td><td>{{.Perc}}%</td><td>{{.FilesWith}}</td><td>{{.FilesWithout}}</td>
			{{end}}

			</tr>	
			{{end}}		
			</tbody>
			</table>
			</table>
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>
			
			
`
	t := template.New("t")

	t, err := t.Parse(htmlTemplate)

	if err != nil {
		panic(err)
	}

	_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

	http.HandleFunc(
		"/", func(w http.ResponseWriter, r *http.Request) {

			err = t.Execute(w, stat)
			if err != nil {
				panic(err)
			}
			go func() {
				defer os.Exit(0)
			}()
		})

	locPort := fmt.Sprintf(":%v", port)
	_ = http.ListenAndServe(locPort, nil)
}

func GenUnique(list []int) []int {
	uniqueSet := make(map[int]bool, len(list))
	for _, x := range list {
		uniqueSet[x] = true
	}
	result := make([]int, 0, len(uniqueSet))
	for x := range uniqueSet {
		result = append(result, x)
	}
	return result
}

func ExportToExcel(snps []bsatstruct.SnpInfo, file string) {

	xlsx := excelize.NewFile()
	// Create a new sheet.
	index := xlsx.NewSheet("SNPs")
	// Set value of a cell.
	var locuses []bsatstruct.SnpInfo

	// Set active sheet of the workbook.
	xlsx.SetActiveSheet(index)
	// Save xlsx file by the given path.

	for _, locus := range snps {

		if locus.Locus != "" {
			locuses = append(locuses, locus)

		}
	}

	for i, snp := range locuses {

		if i > 0 {
			xlsx.SetCellValue("SNPs", fmt.Sprintf("A%v", i+1), snp.Locus)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("B%v", i+1), snp.Name)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("C%v", i+1), snp.APos)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("D%v", i+1), fmt.Sprintf("%v%v>%v", snp.PosInGene, snp.NucInPosCoding, snp.Alt))
			xlsx.SetCellValue("SNPs", fmt.Sprintf("E%v", i+1), fmt.Sprintf("%v/%v", snp.RefCodon, snp.AltCodon))
			xlsx.SetCellValue("SNPs", fmt.Sprintf("F%v", i+1), fmt.Sprintf("%v%v%v", snp.RefAA, snp.CodonNbrInG, snp.AltAA))
			xlsx.SetCellValue("SNPs", fmt.Sprintf("G%v", i+1), snp.Mutation)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("H%v", i+1), snp.Tang)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("I%v", i+1), snp.TangIdxVal)
			xlsx.SetCellValue("SNPs", fmt.Sprintf("J%v", i+1), snp.Product)

		} else {
			xlsx.SetCellValue("SNPs", "A1", "Locus")
			xlsx.SetCellValue("SNPs", "B1", "Gene")
			xlsx.SetCellValue("SNPs", "C1", "Position")
			xlsx.SetCellValue("SNPs", "D1", "Mutation")
			xlsx.SetCellValue("SNPs", "E1", "Codons")
			xlsx.SetCellValue("SNPs", "F1", "AA")
			xlsx.SetCellValue("SNPs", "G1", "Type")
			xlsx.SetCellValue("SNPs", "H1", "Complex Index")
			xlsx.SetCellValue("SNPs", "I1", "Index Sum")
			xlsx.SetCellValue("SNPs", "J1", "Product")

		}

	}
	err := xlsx.SaveAs(file)
	if err != nil {
		fmt.Println(err)
	} else {
		fmt.Printf("The %v file was created successfully\n", file)
	}
}

func MatrixPrint(headers strings.Builder, buffer strings.Builder, fileOut string) {
	pterm.EnableDebugMessages()

	if buffer.Len() != 0 && headers.Len() != 0 {
		fOut, err := os.Create(fileOut)
		if err != nil {
			log.Fatal("Cannot create file", err)
		}
		defer func(fOut *os.File) {
			err := fOut.Close()
			if err != nil {
				log.Fatal(err)
			}
		}(fOut)
		_, _ = fmt.Fprintf(fOut, headers.String())
		_, _ = fmt.Fprintf(fOut, buffer.String())

		// fmt.Printf("\n\nWell done!\n")
		pterm.Success.Printf("%v saved\n", fileOut)
	}
}

func OpenBrowser(url string) {
	var err error

	switch runtime.GOOS {
	case "linux":
		err = exec.Command("xdg-open", url).Start()
	case "windows":
		err = exec.Command("rundll32", "url.dll,FileProtocolHandler", url).Start()
	case "darwin":
		err = exec.Command("open", url).Start()
	default:
		err = fmt.Errorf("unsupported platform")
	}
	if err != nil {
		log.Fatal(err)
	}

}

func PrintSeqRange(rangeList []bsatstruct.RangePosInfo, web bool, port string, toCircos bool) {
	var basicAnnotation string
	switch web {
	case false:
		if toCircos == true && bsatstruct.Flag.StatCircosGC == false {
			basicAnnotation = "{{range $element := .}}" +
				"{{.Genome}}\t{{.Start}}\t{{.End}}\t{{.GeneName}}\n" +
				"{{end}}"

		} else if toCircos == true && bsatstruct.Flag.StatCircosGC == true {
			basicAnnotation = "{{range $element := .}}" +
				"{{if (ge .GC 1.0) }}" +
				"{{.Genome}}\t{{.Start}}\t{{.End}}\t{{.GC}}\n" +
				"{{end}}" +
				"{{end}}"
		} else if bsatstruct.Flag.StatMakeSeq == true || bsatstruct.Flag.StatMkSeq == true {

			basicAnnotation = "{{range $element := .}}" +
				">" +
				"{{.Genome}}[{{.Start}}" + ":" + "{{.End}}]{{.GeneName}}\n" +
				"{{.Seq}}\n" +
				"{{end}}"

		} else {
			basicAnnotation = "{{range $element := .}}" +
				"{{if .Seq}}" +
				"{{.GeneName}}\t{{.Start}}\t{{.End}}\t{{.Len}}\t{{.GC}}\t{{.Seq}}\t{{.Prod}}\t{{.Doubles}}\n" +
				"{{else}}" +
				"{{.GeneName}}\t{{.Start}}\t{{.End}}\t{{.Len}}\t{{.GC}}\t{{.Prod}}\t{{.Doubles}}\n" +
				"{{end}}" +
				"{{end}}"

		}

		t := template.New("basic")
		t, err := t.Parse(basicAnnotation)
		err = t.Execute(os.Stdout, rangeList)
		if err != nil {
			log.Fatal("Parse: ", err)
			return
		}
	case true:

		var htmlAnnotation = `   <!DOCTYPE html>
			<html>
			<head>
			<meta charset="utf-8">			
			</head>		
			<table width="100%" cellspacing="0" cellpadding="4" border="1">
			<body>
			<tr>			
			{{range $element := .}}	
			{{if .Seq}}
			<td>{{.GeneName}}</td><td>{{.Start}}</td><td>{{.End}}</td><td>{{.Len}}</td><td><p title="GC content: {{.GC}}"><textarea rows="3" style="width:400px; word-wrap:break-word;">{{.Seq}}</textarea></p></td><td>{{.Len}}</td><td><p title="{{.Note}}">{{.Prod}}</p></td><td>{{.Doubles}}</td>
			{{else}} 
			<td>{{.GeneName}}</td><td>{{.Start}}</td><td>{{.End}}</td><td>{{.Len}}</td><td><p title="{{.Note}}">{{.Prod}}</p></td><td>{{.Doubles}}</td>
			{{end}}
			</tr>	
			{{end}}		
			</body>
			</table>
		
			<table width="100%" cellspacing="0" cellpadding="4" border="0">
			<tr>
			<td><a href="http://bsatool.ru" target="_blank">Created by BSATool (Bacterial Snp Annotation Tool)</a></td>
			</tr>
			</table>
			
			
`
		tH := template.New("html")

		tH, err := tH.Parse(htmlAnnotation)

		if err != nil {
			panic(err)
		}

		_ = browser.OpenURL(fmt.Sprintf("localhost:%v/", port))

		http.HandleFunc(
			"/", func(w http.ResponseWriter, r *http.Request) {

				err = tH.Execute(w, rangeList)
				if err != nil {
					panic(err)
				}
				go func() {
					defer os.Exit(0)
				}()
			})

		locPort := fmt.Sprintf(":%v", port)
		_ = http.ListenAndServe(locPort, nil)

	}
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

func ChkPosExists(s []bsatstruct.AllPositionsInGene, pos int, alt string) bool {
	for _, a := range s {
		if a.Pos == pos && a.Alt == alt {
			return true
		}
	}
	return false
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

func GetGenePosByName(locus string) (int, int) {
	var start, end int

	g := bsatstruct.GenePositions[locus]

	if g.Start != 0 {
		start = g.Start
		end = g.End
	}

	return start, end
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

func GetListOfGenes() []string {
	var (
		allGenes []string
	)

	for _, g := range bsatstruct.FullGenesInfo {
		allGenes = append(allGenes, g.Locus)

	}

	return allGenes
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

// func GetCoordRange(start, end int) {
//
// 	var coordArray []int
//
// 	coordArray = append(coordArray, start)
//
// 	for i := start; i <= end; i++ {
//
// 		for _, g := range bsatstruct.FullGenesInfo {
//
// 			if g.Start == i {
//
// 				coordArray = append(coordArray, g.Start)
//
// 			} else if g.End == i {
//
// 				coordArray = append(coordArray, g.End)
//
// 			}
//
// 		}
// 	}
// 	coordArray = append(coordArray, end)
// 	sort.Ints(coordArray)
//
// 	var res []bsatstruct.RangePosInfo
//
// 	for _, val := range coordArray {
//
// 		res = append(res, bsatseq.GetSeqRange(fmt.Sprintf("%v:%v", val, val), bsatstruct.Flag.GbNoSeq))
//
// 	}
//
// 	var last string
// 	for _, val := range res {
// 		if val.GeneName != last {
// 			fmt.Println(val.GeneName, val.Prod)
// 		}
// 		last = val.GeneName
//
// 	}
//
// }

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

func GetRevComplement(sequence string) string {

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

func GetComplement(sequence string) string {

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
