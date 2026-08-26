(* ::Package:: *)

<<Anyonica`
(*testcase*)
fc = FCL[[33]]
fs = FSymbols[fc]
rs = RSymbols[fc]
ps = PSymbols[fc]


ClearAll[minpolycoeffsstring]
minpolycoeffsstring[n_] := minpolycoeffsstring[n] = ToString@CoefficientList[ MinimalPolynomial[n] @ x , x]
floatstring[n_] := ToString[DecimalForm@InfN[n,32]]//StringReplace["I" -> " * im"]
tojuliabrackets[s_String]:=StringReplace[s,{"{"->"[","}"->"]"}];
ruletodictstring[s_[a__]->n_] := 
	Module[ 
		{keystring, tostring,valstring},
		keystring = ToString[{a}] // tojuliabrackets;
		tostring = " => ";
		valstring = "(" <> minpolycoeffsstring[n] <>", " <> floatstring[n] <> ")" // tojuliabrackets;
		StringJoin[ keystring, tostring, valstring ]
	];
	
tojuliadict[ symbols_ ] :=
	StringJoin[
		"Dict(\n",
		StringRiffle[ ruletodictstring /@ symbols, ",\n"],
		"\n)"
	]
	
allsymbolsdict[cat_] :=
	Module[ 
		{ fs, rs, ps, bs },
		fs = FSymbols @ cat; rs = RSymbols @ cat; ps = PSymbols @ cat;
		bs = If[ 
			BraidedQ[cat],
			tojuliadict[rs],
			"Dict()"
		];
		StringJoin[
		"Dict(\n",
		"\"fsymbols\" => ", tojuliadict[fs],",\n",
		"\"rsymbols\" => ", bs, ",\n",
		"\"psymbols\" => ", tojuliadict[ps],"\n",
		")"
		]
	]


expdir = "/home/gert/Projects/Lyctr.jl/src/data/mathematica/Symbols/";

codeTofn[i_]:= StringRiffle[ ToString/@FC[FCL[[i]]], "_" ];

exportsymbols[i_Integer]:=
	Module[{fn},
		fn = FileNameJoin[{expdir,"symbols_"<>codeTofn[i]<>".jl" }];
		Export[fn,allsymbolsdict[FCL[[i]]],"Text"]
	]

exportfs[i_Integer]:= Module[{fn},
	fn = FileNameJoin[{expdir,"fsymbols","fs_"<>ToString[i]<>".jl" }];
	Export[fn, tojuliadict[FSymbols[FCL[[i]]]], "Text" ]
]

exportrs[i_Integer]:= Module[{fn},
	fn = FileNameJoin[{expdir,"rsymbols","rs_"<>ToString[i]<>".jl" }];
	If[ 
		!BraidedQ[FCL[[i]]], 
		Export[fn, tojuliadict[RSymbols[FCL[[i]]]], "Text" ]
	]
]

exportps[i_Integer]:= Module[{fn},
	fn = FileNameJoin[{expdir,"psymbols","ps_"<>ToString[i]<>".jl" }];
	Export[fn, tojuliadict[PSymbols[FCL[[i]]]], "Text" ]
]


PMap[ exportsymbols, Range[889:Length[FCL]] ]


(* ::Section:: *)
(*Exporting Gauge Split Bases*)


(* We will export these using sparse rows which Oscar can use to create the matrices *)
gsts = Import["/home/gert/Projects/Lyctr.jl/src/data/mathematica/GaugeSplitStuff/gauge_split_transforms.mx"];


(* First we check whether these transforms give a gsb *)
applyGST[{sa_,n_}][symbols_]:=
	With[{V = Normal[sa]},
	Map[
      PowerDot[ symbols, Transpose[#] ]&,
      {
        V[[;;,;;n]],
        V[[;;,n+1;;]]
      }
    ]
   ]
   
isGaugeInvariant[ring_][l_List]:= 
	Module[{gtl,g},
		gtl = l/.{\[ScriptCapitalF][i__]:>GaugeTransform[g][\[ScriptCapitalF][i]],\[ScriptCapitalR][i__]:>GaugeTransform[g][\[ScriptCapitalR][i]],\[ScriptP][i__]:>GaugeTransform[ring,g][\[ScriptP][i]]};
		gtl == l
	];
	
checkGaugeInvariance[ i_ ] :=
Module[
	{fc = FCL[[i]],gst = gsts[[i]],r,rsymb,symb},
	r = FusionRing[fc];
	rsymb = If[ BraidedQ[fc], RSymbols[r],{}];
	symb = Join[FSymbols[r],rsymb,PSymbols[r]];
	isGaugeInvariant[r][applyGST[gst][symb][[1]]]
] 


(*It seems we assumed that all transforms keep the vacuum fixed. I don't like this assumption so I'm not going to use these transforms for GSBs in Julia. They can still be used to find smaller fields for 
the F-symbols which is nice*)


checkGaugeInvariance[2]


toOscar[gst_] :=
Module[{arr,ar,rowind,colindval},
	arr = gst[[1]];
	ar = Most@ArrayRules @ arr;
	rowind[{a_,b_}->c_] := a;
	colindval[{a_,b_}->c_]:={b,c};
	<|
		"dimension" -> First @ Dimensions[arr],
		"sparse_rows" -> Map[ colindval, GatherBy[ar,rowind], {2}],
		"n_invariants" -> gst[[2]]
	|>
];


CreatePacletArchive["~/Projects/Anyonica","~/Documents/Paclets/"]



(* ::Section:: *)
(*Exporting modular data and quantum dimensions*)


toqqbid[x_]:= 
	Module[{mp,rn,maxn,y},
		mp = MinimalPolynomial @ x;
		maxn =Exponent[mp@y,y];
		rn =FirstCase[ Range[1,maxn], i_/;RootReduce[x-Root[mp,i]]==0 ];
		StringRiffle[ ToString/@CoefficientList[mp@y,y],
			"_"
		]<>"__"<>ToString[rn]
	]


codestring[cat_] := StringRiffle[ToString/@FC[cat],"_"] 


(* ::Subsection:: *)
(*Modular data*)


smatdict = 
	Association @ 
	PMap[
		codestring[#]-> 
		If[
			ModularQ[#],
			Map[toqqbid,SMatrix[#],{2}],
			Null 
		]&,
		FCL
	];
Export["/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/s_matrices.json",smatdict]


twistdict =  
	Association @ 
	Monitor[Table[
		codestring[cat]-> 
		If[
			ModularQ[cat],
			Map[toqqbid,Simplify[Arg[Values[Twists[cat]]]/(2 Pi)]],
			Null 
		],
		{cat,FCL}
	],cat];


Export["/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/spins.json",twistdict]


(* ::Subsection:: *)
(*Quantum dimensions*)


qddict =  
	Association @ 
	Monitor[Table[
		codestring[cat]-> 
		toqqbid/@Values@QuantumDimensions@cat,
		{cat,FCL}
	],cat];


Export["/home/gert/Projects/Lyctr.jl/src/data/fusion_categories/quantum_dims.json",qddict]
