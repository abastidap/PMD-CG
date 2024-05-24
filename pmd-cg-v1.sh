#!/bin/bash

PROGRAMAS=/home/adolfo/adolfo/investigacion/articulos/idr_ensembles/github/programs

# $1 -> Secuencia
# $2 -> numero de primera estructura
# $3 -> numero de la ultima estructura
# $4 -> resolucion mapas ramachandran
# $5 -> directorio de los tripeptidos
# $6 -> numero de trayectorias para la estadistica de los tripeptidos
# $7 -> ca-ha (calculo de solapamientos por distancias), chimera
# $8 -> scwrl4 o faspr o none (en este caso se usa la estructura generada por pymol)
# $9 -> rotar los angulos chi1 (chiyes/chino)
# $10 -> tipo de probabilidad p1, p2izq o p2de

seq=AHSSHLKSKKGQSTSRHKKL
ntrayi=$1
ntrayf=$2
res=360
dirt=$3
ntri=100
fil=chimera
side=scwrl4
chi=chino
prob=p1
PROGRAMAS=$4
version=1

cat <<EOF >AHSSHLKSKKGQSTSRHKKL.in
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
betap
-59.,169.,50.,50.
betapr
60.,168.,50.,50.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
2
alpha
-90.,-35.,90.,85.
gamma
90.,0.,90.,180.
EOF

((ntray=$ntrayf-$ntrayi+1))

local=$(pwd)

outputdir=$local/pmd-cg-v$version/$seq
mkdir -p $local/pmd-cg-v$version/$seq

rm -fr borrar
mkdir -p borrar

padding=5
zi=""
for ((a=1; a <= "$padding" ; a++))
do
zi=$zi"0"
done


# almacenamos cada residuo de la cadena en el vector r

len=$(echo $seq |awk '{print length}')

for ((a=1; a <= "$len" ; a++))
do
    r[$a]=$(echo $seq | head -c $a | tail -c 1)
    # echo "Res "$a" = "${r[$a]}
done

# comprobamos si existen los calculos de los tripeptidos y copiamos los ficheros raman, de distancias y los de omega

rm -rf borrar
mkdir -p borrar

((len1=$len-1))
((len2=$len-2))


for ((a=1; a <= "$len2" ; a++))
do
    ((b=$a+1))
    ((c=$a+2))

    t=${r[$a]}${r[$b]}${r[$c]}
    
    dir=$dirt/$t/test_equilibrado/charmm36/geometria_extendida_xtc/datos-10ns-exp10fs-eq1+1ns/distancia_esqueleto-v1-tray-1-$ntri


    if [ -d "$dir" ]; then
	# echo "si existe directorio de distancias de "$t
	cp $dir/histo-$t-*.dat borrar
	else
	echo "no existe directorio distancias de "$t
	exit
    fi

    dir=$dirt/$t/test_equilibrado/charmm36/geometria_extendida_xtc/datos-10ns-exp10fs-eq1+1ns/rama-v2-tray-1-$ntri/corr_cortes/res-$res


    if [ -d "$dir" ]; then
	# echo "si existe directorio raman de "$t
	cp $dir/raman-$t-di1*.dat borrar
	cp $dir/raman-$t-di2*.dat borrar
	cp $dir/raman-$t-di3*.dat borrar
	else
	echo "no existe directorio raman de "$t
	exit
    fi
    
    dir=$dirt/$t/test_equilibrado/charmm36/geometria_extendida_xtc/datos-10ns-exp10fs-eq1+1ns/diedros-v4-tray-1-$ntri


    if [ -d "$dir" ]; then
	# echo "si existe directorio omega de "$t
	cp $dir/histo-omega-res-1-2.dat borrar/histo-omega-res-1-2-$t.dat
	cp $dir/histo-omega-res-2-3.dat borrar/histo-omega-res-2-3-$t.dat
    else
	echo "no existe directorio omega de "$t
	exit
    fi


    if [ -d "$dir" ]; then

	if [ "$a" == "1" ];then
	    if [ "${r[$a]}" != "A" ] && [ "${r[$a]}" != "G" ] && [ "${r[$a]}" != "P" ]; then
		cp $dir/histo-chi1-res-1.dat borrar/histo-chi1-res-1.dat
	    fi
	fi

	if [ "${r[$b]}" != "A" ] && [ "${r[$b]}" != "G" ] && [ "${r[$b]}" != "P" ]; then
	    cp $dir/histo-chi1-res-2.dat borrar/histo-chi1-res-$b.dat
	fi

	if [ "$a" == "$len2" ];then
	    if [ "${r[$a]}" != "A" ] && [ "${r[$a]}" != "G" ] && [ "${r[$a]}" != "P" ]; then
		cp $dir/histo-chi1-res-3.dat borrar/histo-chi1-res-$len.dat
	    fi
	fi
	
    else
	echo "no existe directorio chi de "$t
	exit
    fi

    
done

# # comprobar si estan todos

# ls -lita borrar/histo-chi1-res*.dat

# # exit 0

cd borrar

# generamos los ficheros con los angulos diedros aleatorios
# generamos un multiplo de angulos porque descartaremos un % de estructuras con solapamientos

((ntraym=$ntray*6))

cat <<EOF >input.in
$seq
$res
$ntraym
EOF

cat ../$seq.in >>input.in

$PROGRAMAS/probabilidades_ramachandran/probabilidades_ramachandran.exe < input.in

mv fort.40 diedros-estadistica.dat
mv fort.50 diedros1.dat
mv fort.51 diedros2.dat
mv fort.52 diedros3.dat

# generamos los ficheros con los valores de omega aleatorios

for ((a=1; a <= "$len2" ; a++))
do
    ((b=$a+1))
    ((c=$a+2))

    t=${r[$a]}${r[$b]}${r[$c]}

    nf=$(cat histo-omega-res-1-2-$t.dat | wc -l)
    ((nf=$nf-1))

    
    cat <<EOF >input.in
$nf
histo-omega-res-1-2-$t.dat
$ntraym
EOF
    $PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
    mv fort.40 omega-$a.dat

done

a=$len2
((b=$a+1))
((c=$a+2))

t=${r[$a]}${r[$b]}${r[$c]}

nf=$(cat histo-omega-res-2-3-$t.dat | wc -l)
    ((nf=$nf-1))
    
cat <<EOF >input.in
$nf
histo-omega-res-2-3-$t.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 omega-$len1.dat

# generamos los ficheros con las distancias aleatorias

t=${r[1]}${r[2]}${r[3]}

nf=$(cat histo-$t-CH31-C1.dat | wc -l)

cat <<EOF >input.in
$nf
histo-$t-CH31-C1.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-CH31-C1.dat

cat <<EOF >input.in
$nf
histo-$t-C1-N1.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-C1-N1.dat

cat <<EOF >input.in
$nf
histo-$t-N1-CA1.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-N1-CA1.dat

cat <<EOF >input.in
$nf
histo-$t-CA1-C2.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-CA1-C2.dat


for ((a=2; a <= "$len1" ; a++))
do
    ((m=$a-1))
    ((p=$a+1))
    t=${r[$m]}${r[$a]}${r[$p]}

    cat <<EOF >input.in
$nf
histo-$t-C2-N2.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-C$a-N$a.dat
    
    cat <<EOF >input.in
$nf
histo-$t-N2-CA2.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-N$a-CA$a.dat
    
    cat <<EOF >input.in
$nf
histo-$t-CA2-C3.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-CA$a-C$p.dat    
    
done

((m=$len-2))
((a=$len-1))
p=$len
((q=$len+1))

t=${r[m]}${r[a]}${r[p]}

    cat <<EOF >input.in
$nf
histo-$t-C3-N3.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-C$p-N$p.dat
    
    cat <<EOF >input.in
$nf
histo-$t-N3-CA3.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-N$p-CA$p.dat
    
    cat <<EOF >input.in
$nf
histo-$t-CA3-C4.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-CA$p-C$q.dat    

cat <<EOF >input.in
$nf
histo-$t-C4-N4.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-C$q-N$q.dat

cat <<EOF >input.in
$nf
histo-$t-N4-CH32.dat
$ntraym
EOF
$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
mv fort.40 dis-N$q-CH32.dat

# generamos los ficheros con los valores de chi aleatorios

for ((a=1; a <= "$len" ; a++))
do
    if [ "${r[$a]}" != "A" ] && [ "${r[$a]}" != "G" ] && [ "${r[$a]}" != "P" ]; then
	nf=$(cat histo-chi1-res-$a.dat | wc -l)
    
	cat <<EOF >input.in
$nf
histo-chi1-res-$a.dat
$ntraym
EOF
	$PROGRAMAS/probabilidades_distancias/probabilidades_distancias.exe < input.in
	mv fort.40 chi1-$a.dat
    fi

done


timei=$(date)
    
if [ "$prob" == "p1" ]; then
    e="1"
fi

if [ "$prob" == "p2izq" ]; then
    e="2"
fi

if [ "$prob" == "p2der" ]; then
    e="3"
fi

    
nval=0
    
# rm -f angulos-aleatorios.dat angulos-reales.dat
    
# touch angulos-aleatorios.dat
# touch angulos-reales.dat
    
nsol=0
    
# creamos el fichero pdb con pymol

# /home/adolfo/adolfo/investigacion/liquidos/construir_seq/v2/construir_seq-v2.sh $seq

cat <<EOF >aux.pml
from pymol import cmd, editor
cmd.set('retain_order', 0)
cmd.fab('$seq', ss=3)
editor.attach_amino_acid("last name C", 'nme')
editor.attach_amino_acid("first name N", 'ace')
save $seq.pdb
EOF

pymol -c aux.pml &> pymol.log

# substituciones para charmm36m ff

sed -i '/TER/d' $seq.pdb


sed -i 's/CT3/NME/' $seq.pdb

sed -i 's/3HA  GLY/1HA  GLY/' $seq.pdb
sed -i 's/ HA  GLY/2HA  GLY/' $seq.pdb

sed -i 's/2HB  LEU/HB1  LEU/' $seq.pdb
sed -i 's/3HB  LEU/HB2  LEU/' $seq.pdb

sed -i 's/2HB  GLN/1HB  GLN/' $seq.pdb
sed -i 's/3HB  GLN/2HB  GLN/' $seq.pdb
sed -i 's/2HG  GLN/1HG  GLN/' $seq.pdb
sed -i 's/3HG  GLN/2HG  GLN/' $seq.pdb

sed -i 's/2HB  ARG/1HB  ARG/' $seq.pdb
sed -i 's/3HB  ARG/2HB  ARG/' $seq.pdb
sed -i 's/2HG  ARG/1HG  ARG/' $seq.pdb
sed -i 's/3HG  ARG/2HG  ARG/' $seq.pdb
sed -i 's/2HD  ARG/1HD  ARG/' $seq.pdb
sed -i 's/3HD  ARG/2HD  ARG/' $seq.pdb

sed -i 's/2HB  ASN/HB1  ASN/' $seq.pdb
sed -i 's/3HB  ASN/HB2  ASN/' $seq.pdb

sed -i 's/2HB  ASP/1HB  ASP/' $seq.pdb
sed -i 's/3HB  ASP/2HB  ASP/' $seq.pdb

sed -i 's/2HB  CYS/1HB  CYS/' $seq.pdb
sed -i 's/3HB  CYS/2HB  CYS/' $seq.pdb
sed -i 's/ HG  CYS/1HG  CYS/' $seq.pdb

sed -i 's/2HB  GLU/1HB  GLU/' $seq.pdb
sed -i 's/3HB  GLU/2HB  GLU/' $seq.pdb
sed -i 's/2HG  GLU/1HG  GLU/' $seq.pdb
sed -i 's/3HG  GLU/2HG  GLU/' $seq.pdb

sed -i 's/2HB  HIS/1HB  HIS/' $seq.pdb    
sed -i 's/3HB  HIS/2HB  HIS/' $seq.pdb    

sed -i 's/2HG1 ILE/1HG1 ILE/' $seq.pdb    
sed -i 's/3HG1 ILE/2HG1 ILE/' $seq.pdb    

sed -i 's/2HB  LYS/1HB  LYS/' $seq.pdb
sed -i 's/3HB  LYS/2HB  LYS/' $seq.pdb
sed -i 's/2HG  LYS/1HG  LYS/' $seq.pdb
sed -i 's/3HG  LYS/2HG  LYS/' $seq.pdb
sed -i 's/2HD  LYS/1HD  LYS/' $seq.pdb
sed -i 's/3HD  LYS/2HD  LYS/' $seq.pdb
sed -i 's/2HE  LYS/1HE  LYS/' $seq.pdb
sed -i 's/3HE  LYS/2HE  LYS/' $seq.pdb

sed -i 's/2HB  MET/1HB  MET/' $seq.pdb
sed -i 's/3HB  MET/2HB  MET/' $seq.pdb
sed -i 's/2HG  MET/1HG  MET/' $seq.pdb
sed -i 's/3HG  MET/2HG  MET/' $seq.pdb

sed -i 's/2HB  PHE/1HB  PHE/' $seq.pdb
sed -i 's/3HB  PHE/2HB  PHE/' $seq.pdb

sed -i 's/2HB  PRO/1HB  PRO/' $seq.pdb
sed -i 's/3HB  PRO/2HB  PRO/' $seq.pdb
sed -i 's/2HG  PRO/1HG  PRO/' $seq.pdb
sed -i 's/3HG  PRO/2HG  PRO/' $seq.pdb
sed -i 's/2HD  PRO/1HD  PRO/' $seq.pdb
sed -i 's/3HD  PRO/2HD  PRO/' $seq.pdb

sed -i 's/2HB  SER/1HB  SER/' $seq.pdb
sed -i 's/3HB  SER/2HB  SER/' $seq.pdb
sed -i 's/ HG  SER/1HG  SER/' $seq.pdb

sed -i 's/2HB  TRP/1HB  TRP/' $seq.pdb
sed -i 's/3HB  TRP/2HB  TRP/' $seq.pdb

sed -i 's/2HB  TYR/1HB  TYR/' $seq.pdb
sed -i 's/3HB  TYR/2HB  TYR/' $seq.pdb

cp $seq.pdb $seq.seg

# creamos el fichero gro

echo "3 4" | gmx pdb2gmx -f $seq.pdb -o $seq.gro -p $seq.top -ter -ff charmm36-mar2019 -water tip3p -ignh &> pdb2gmx.log

rm -f *posre.itp* aux.pml

# contamos el numero de atomos como control porque scwrl4 parece añadir en casos especiales un atomo adicional para HIS

nc=$(sed '2q;d' $seq.gro)

# echo "Numero de atomos "$nc
    
# generar topologia

gmx editconf -f $seq.gro -o $seq-box.gro -bt cubic -box 7.0 7.0 7.0 &> /dev/null 

echo "3 4" | gmx pdb2gmx -f $seq-box.gro -o $seq-boxb.gro -p $seq-boxb.top -ter -ff charmm36-mar2019 -water tip3p &> /dev/null 



# copiamos un fichero mdp

t=${r[1]}${r[2]}${r[3]}

cp $dirt/$t/equilibrado/charmm36/geometria_extendida_xtc/equiNVT.mdp .

# cp equiNVT.mdp $local/generacion_estructuras-v$version
    
# creamos un .tpr (warning si la carga total no es nula)

gmx grompp -f equiNVT.mdp -c $seq-boxb.gro -p $seq-boxb.top -o $seq-boxb.tpr -maxwarn 1 &> /dev/null 

    
# bucle de conformaciones a generar

b=0
   
while [ "$nval" -lt "$ntray" ]
do

    cp $seq.seg $seq.pdb
    
    ((b=$b+1))
 
    ((i=$len*(b-1)+1))
    ((f=$i+$len-1))

    sed -n "$i,$f p"  diedros$e.dat > dihedrals.in 

    # /home/adolfo/adolfo/investigacion/liquidos/construir_seq/v2/construir_seq-v2.sh $seq
    
    # modificamos los valores de phi, psi y omega con pymol
		
    echo "load "$seq".pdb" >  aux.pml

    for ((a=1; a <= "$len" ; a++))
    do

    	((r1=$a-1))
    	r2=$a
    	((r3=$a+1))

    	phi=$(awk {'print $1'}  dihedrals.in | head -$a | tail -1)
    	psi=$(awk {'print $2'}  dihedrals.in | head -$a | tail -1)
	
    	echo "set_dihedral resi "$r1" and name C, resi "$r2" and name N, resi "$r2" and name CA, resi "$r2" and name C, "$phi >> aux.pml
    	echo "set_dihedral resi "$r2" and name N, resi "$r2" and name CA, resi "$r2" and name C, resi "$r3" and name N, "$psi >> aux.pml
	
    done

    
    for ((a=1; a <= "$len1" ; a++))
    do
	((r1=$a))
	((r2=$a+1))
	    
	ang=$(sed "${b}q;d" omega-$a.dat)
	echo "set_dihedral resi "$r1" and name CA, resi "$r1" and name C, resi "$r2" and name N, resi "$r2" and name CA, "$ang >> aux.pml
		    
    done
		
    echo "save aux.pdb" >> aux.pml
    pymol -c aux.pml &> pymol.log

    mv aux.pdb $seq.pdb

    sed -i '/TER/d' $seq.pdb
   
    # modificamos las distancias del esqueleto con chimera
		
    cat <<EOF > input.com
open $seq.pdb
EOF
		
    dis=$(sed "${b}q;d" dis-CH31-C1.dat)
    dis=$(echo "scale=3; $dis * 10" | bc)

    echo "adjust length $dis :0@CH3 :0@C" >> input.com

    ((lenm1=$len-1))
    ((lenp1=$len+1))
    for ((a=1; a <= "$lenp1" ; a++))
    do
	dis=$(sed "${b}q;d" dis-C$a-N$a.dat)
	dis=$(echo "scale=3; $dis * 10" | bc)
	((m=$a-1))
	echo "adjust length $dis :$m@C :$a@N" >> input.com
    done

    for ((a=1; a <= "$len" ; a++))
    do
	dis=$(sed "${b}q;d" dis-N$a-CA$a.dat)
	dis=$(echo "scale=3; $dis * 10" | bc)
	echo "adjust length $dis :$a@N :$a@CA" >> input.com
    done

    for ((a=1; a <= "$len" ; a++))
    do
	((m=$a+1))
		    
	dis=$(sed "${b}q;d" dis-CA$a-C$m.dat)
	dis=$(echo "scale=3; $dis * 10" | bc)
	echo "adjust length $dis :$a@CA :$a@C" >> input.com
		    
    done

    dis=$(sed "${b}q;d" dis-N$lenp1-CH32.dat)
    dis=$(echo "scale=3; $dis * 10" | bc)
    # echo "adjust length $dis :$lenp1.het@N :$lenp1.het@CH3" >> input.com
    echo "adjust length $dis :$lenp1@N :$lenp1@CH3" >> input.com

    echo "write format pdb #0 $seq.pdb" >> input.com

    chimera --nogui input.com &> chimera.log

    sed -i 's/HETATM/ATOM  /' $seq.pdb
    sed -i '/CONECT/d' $seq.pdb 
    sed -i '/END/d' $seq.pdb 

    if [ "$side" != "none" ]; then
    
	# store terminal groups which will be deleted by scwrl4

	mv $seq.pdb aux.pdb

	grep ACE aux.pdb > head.pdb
	grep NME aux.pdb > tail.pdb

	# extract backbone atoms

	rm -f aux2.pdb
	
	l=$(wc -l < aux.pdb)
    
	for ((a=1; a <= "$l" ; a++))
	do
   
	    c=$(head -$a aux.pdb | tail -1 | awk {'print $1'})
    
	    if [ $c == "ATOM" ]; then

		d=$(head -$a aux.pdb | tail -1 | awk {'print $4'})

		if [ $d != "ACE" ]; then
		    
		    if [ $d != "NME" ]; then

			sed "${a}q;d" aux.pdb >> aux2.pdb
		
		    fi
		fi
		
	    fi
	done

	# use scwrl4 or faspr to add lateral chains

	if [ "$side" == "scwrl4" ]; then   
	    $SCWRL4/Scwrl4 -t -i aux.pdb -o aux2.pdb &> scwrl4.log
	fi

	if [ "$side" == "faspr" ]; then		
	    sed -i '/ACE/d' aux.pdb
	    sed -i '/NME/d' aux.pdb
	    sed -i '/TER/d' aux.pdb
	    $FASPR/FASPR -i aux.pdb -o aux2.pdb &> faspr.log
	    sed -i -e "1d" aux2.pdb		
	fi
	    
	cat head.pdb aux2.pdb tail.pdb > $seq.pdb

	sed -i '/TER/d' $seq.pdb

    fi

    # # comprobar la numeracion de residuos en el pdb

    # cat $seq.pdb
    
   
# modificamos los angulos chi1 con pymol

    if [ "$chi" == "chiyes" ]; then   
		
	echo "load "$seq".pdb" >  aux.pml

	for ((a=1; a <= "$len" ; a++))
	do
	    if [ "${r[$a]}" != "A" ] && [ "${r[$a]}" != "G" ] && [ "${r[$a]}" != "P" ]; then
	
		ang=$(sed "${b}q;d" chi1-$a.dat)

		nom4="CG"
		if  [ "${r[$a]}" == "C" ]; then
		    nom4="SG"
		fi
		if  [ "${r[$a]}" == "I" ] ||  [ "${r[$a]}" == "V" ] ; then
		    nom4="CG1"
		fi
		if  [ "${r[$a]}" == "S" ]; then
		    nom4="OG"
		fi
		if  [ "${r[$a]}" == "T" ]; then
		    nom4="OG1"
		fi
		
		echo "set_dihedral resi "$a" and name N, resi "$a" and name CA, resi "$a" and name CB, resi "$a" and name "$nom4", "$ang >> aux.pml

	    fi	
	done
		
	echo "save aux.pdb" >> aux.pml
	pymol -c aux.pml &> pymol.log

	# # comprobar rotacion chi
	# cat aux.pml
	# exit 0
	
	mv aux.pdb $seq.pdb

	sed -i '/TER/d' $seq.pdb

	
    fi

# fase final
    
    echo "3 4" | gmx pdb2gmx -f $seq.pdb -o $seq.gro -p borrar.top -ter -ff charmm36-mar2019 -water tip3p -ignh &> pdb2gmx.log

    rm *posre.itp* *borrar.top*

    sed -i '$ d' $seq.gro
    echo "7.0000 7.0000 7.0000 " >> $seq.gro

    # busqueda de solapamientos (clashes)
	
    if [ $fil == "ca-ha" ]; then
	
	# comprobamos si hay solapamiento entre carbonos alpha no consecutivos

	tail -n +2 $seq.gro > aux.dat

	grep '\<CA\>' aux.dat |awk '{print $4" "$5" "$6}' > prueba.dat

	n=$(wc -l < prueba.dat)

	cat <<EOF >aux.in
$n	
prueba.dat
0.4		
EOF
	$PROGRAMAS/distancias_ca_filtro/distancias_ca_filtro.exe < aux.in
	
	m=$(awk '{print $1}' fort.40)

	t=$(awk '{print $1" "$2" "$3}' fort.40)
	
	cat <<EOF >aux.in
$seq.gro
0.2		
EOF
	
	$PROGRAMAS/distancias_heavy_atoms/distancias_heavy_atoms.exe < aux.in

	n=$(awk '{print $1}' fort.40)

	u=$(awk '{print $1" "$2" "$3}' fort.40)

	if [ $m == "No" ] && [ $n == "No" ]; then

	    clash="No"
	    	    
	else
		
	    clash="Si"
	    echo "Solapamiento ca "$t
	    echo "Solapamiento ha "$u
	    ((nsol=$nsol+1))
	    echo "(total err "$nsol")"
		
	fi
	    
    fi

    if [ $fil == "chimera" ]; then

	cat <<EOF >aux.com
open $seq.pdb
findclash #0 test self intraRes true saveFile contacts.log
EOF

	chimera --nogui aux.com &> chimera.log

	con=$(sed "7q;d" contacts.log | awk '{print $1}')
	    
	if [ "$con" == "0" ]; then

	    clash="No"

	else

	    clash="Si"
	    ((nsol=$nsol+1))
	    echo "Overlaps chimera "$con" (total err "$nsol")"
	
	fi
	    
    fi

	
# comprobamos si el numero de atomos es correcto

    na=$(sed '2q;d' $seq.gro)

    if [ "$na" == "$nc" ];then

	if [ "$clash" == "No" ]; then
	    
	    # gmx rama -f $seq.gro -s $seq-boxb.tpr -xvg none &> /dev/null

	    # cat dihedrals.in >> angulos-aleatorios.dat
	    # awk '{print $1" "$2}' rama.xvg >> angulos-reales.dat
	    # rm rama.xvg   
	    
	    ((nval=$nval+1))

	    ((num=$ntrayi+$nval-1))
	    la=${#num}
	    let "cla=$padding-$la"
	    z=$zi
	    c=${z:$la:$cla}$num
		
	    mkdir -p $outputdir/tray${c:0:1}
	    mkdir -p $outputdir/tray${c:0:1}/tray$c
		
	    cp $seq.gro $outputdir/tray${c:0:1}/tray$c/$seq-eqNVEe$c.gro
	    gmx convert-trj -f $outputdir/tray${c:0:1}/tray$c/$seq-eqNVEe$c.gro -o $outputdir/tray${c:0:1}/tray$c/$seq-eqNVEe$c.xtc &> /dev/null
	    gmx editconf -f $outputdir/tray${c:0:1}/tray$c/$seq-eqNVEe$c.gro -o $outputdir/tray${c:0:1}/tray$c/$seq-eqNVEe$c.pdb &> /dev/null

    
	    echo "Structure "$num" ("$ntrayi-$ntrayf")" 

# con $prob"-"$fil"-"$side"-"$chi

	else
		
	    rm *$seq.gro* *$seq.pdb*
	    
	fi
	    
    else
	echo "Numero atomos correcto "$nc" actual "$na
	rm *$seq.gro* *$seq.pdb*
	    
    fi

	
done


echo " " >> $outputdir/info-$ntrayi-$ntrayf.dat
echo "Estructuras con solapamiento = "$nsol >> $outputdir/info-$ntrayi-$ntrayf.dat

    
# # calculamos las correlaciones entre los angulos aleatorios y los generados por el algoritmo

# cat <<EOF >auxr.in
# auxr.dat      
# EOF

  
# echo " " >>  $outputdir/info-$ntrayi-$ntrayf.dat
# echo "Datos de correlacion entre angulos aleatorios y reales " >>  $outputdir/info-$ntrayi-$ntrayf.dat
# echo " " >> $outputdir/info-$ntrayi-$ntrayf.dat


# for ((a=1; a <= "$len" ; a++))
# do

#     sed -n "$a~$len p"  angulos-aleatorios.dat > a.dat
#     sed -n "$a~$len p"  angulos-reales.dat > b.dat
    
    
#     awk '{print $1}' a.dat > a1.dat
#     awk '{print $1}' b.dat > b1.dat

#     paste -d" " a1.dat b1.dat > auxr.dat

#     $PROGRAMAS/regresion_lineal/regresion_lineal.exe < auxr.in

#     cat fort.40 >> $outputdir/info-$ntrayi-$ntrayf.dat 

#     awk '{print $2}' a.dat > a2.dat
#     awk '{print $2}' b.dat > b2.dat

#     paste -d" " a2.dat b2.dat > auxr.dat

#     $PROGRAMAS/regresion_lineal/regresion_lineal.exe < auxr.in

#     cat fort.40 >>  $outputdir/info-$ntrayi-$ntrayf.dat
    
# done

timef=$(date)

echo "========================================================" >> $outputdir/info-$ntrayi-$ntrayf.dat
echo "inicio "$timei >> $outputdir/info-$ntrayi-$ntrayf.dat
echo "fin    "$timef >> $outputdir/info-$ntrayi-$ntrayf.dat
echo "========================================================" >> $outputdir/info-$ntrayi-$ntrayf.dat


cd ..
rm -r borrar
rm AHSSHLKSKKGQSTSRHKKL.in
