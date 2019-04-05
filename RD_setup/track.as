table CNV_TRACK
"Copy number track "
(
string  chrom;          "Reference sequence chromosome or scaffold"
uint    chromStart;     "Start position of feature on chromosome"
uint    chromEnd;       "End position of feature on chromosome"
string  name;           "Individual name"
uint    score;          "Score"
char[1] strand;         "+ or - for strand"
uint    thickStart;     "Coding region start"
uint    thickEnd;       "Coding region end"
uint reserved;	        "LightGray:0, DarkGray:1, Black:2, MediumBlue:3, RoyalBlue:4, Periwinkle:5, LightGreen:6, Yellow:7, Orange:8, Maroon:9, Red:10+"
string  ID;    	        "Copy Number"
)