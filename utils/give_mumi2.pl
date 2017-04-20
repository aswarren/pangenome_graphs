#! /usr/bin/perl -w

use strict;
use Getopt::Long;

=head1  NAME

give_mumi.pl

=head1 DESCRIPTION

A program that calculates the distance between two genomes. 
Entry file: MUMmer3 output of a mummer3 run with options -b -c -l 19
Example:
mummer -mum -b -c -l 19 U00096_GR.fsa AE014075_GR.fsa >MUM_K12_CFT

First it cancels the included MUMs in the two genomes.
Second, it cancels the double_overlapping MUMs in the two genomes
Third, it trims the overlapping MUMs in the genome 0 
Fourth, cancels the included MUMs on genome 1 created by the preceding step
Fifth,it trims the overlapping MUMs in the genome 1 
The last thing is the calculation of the distance with the formula : d(G0,G1)=(Sum(Non_Overlapping_MUMs_Lengths))/(Mean_Genome_Lengths).

Then treats G1 before G0 

Returns to the screen MUMI average
=head1  USAGE

example
./give_mumi.pl MUM_K12_CFT_19 -l1 4639675 -l2 5231428 

Parameters -l1 and -l2 are obligatory, correspond to the length of the two genomes. 
=cut 


#Gestion des paramètres passes en ligne de commande


my ($gen1_len, $gen2_len);
GetOptions( "l1=s" => \$gen1_len,
	    "l2=s" => \$gen2_len);

my $ficName = $ARGV[0];	    


#Parcours fichier Mummer
=head2 Read_mums
  
  Description:
  
    Prend en entrée le nom de fichier mummer
    Retourne 

   
  Arguments:
  
    $ficName 	Prend en entrée le nom de fichier mummer
  Returns:
  
    $rtab_mums    Liste des mums. 
    Pour chacun, liste de deux dictionnaires clés G0,  G1, et ident (pour identifiant du mum)
    Pour les clés G0 et G1, valeur est un dictionnaire ayant pour clés init, fin, long, sens (vaut 1 ou 2)
    
=cut 

sub Read_mums {
    my $ficName=shift;
    open(FIC, $ficName) || die "can't open file $ficName";
    my $i=0;
    my $nb_mums_tot=0;
    my $nb_mums_dir=0;
    my $nb_mums_rev=0;
    my $flag_reverse=0;
    my @tab_mums=();
while (my $ligne=<FIC>){
    chomp($ligne);
    if ($ligne=~/^>/ && $ligne!~/Reverse/){
    	$flag_reverse=0;
    }elsif($ligne=~/^>/ && $ligne=~/Reverse/){
    	$flag_reverse=1;
    }else{      
    
	$nb_mums_tot++;
	
	if(($ligne=~/(\d+)\s+(\d+)\s+(\d+)/) && $flag_reverse==0){
	    $nb_mums_dir++;
	    
	    $tab_mums[$i]->{G0}->{sens}=1;
	    $tab_mums[$i]->{G0}->{init}=$1;
	    $tab_mums[$i]->{G0}->{long}=$3;
	    $tab_mums[$i]->{G0}->{fin}=$1+$3-1;
	    $tab_mums[$i]->{G1}->{sens}=1;
	    $tab_mums[$i]->{G1}->{init}=$2;
	    $tab_mums[$i]->{G1}->{long}=$3;
	    $tab_mums[$i]->{G1}->{fin}=$2+$3-1;	 
	    $tab_mums[$i]->{ident}=$nb_mums_tot;
	    
	}
	elsif(($ligne=~/(\d+)\s+(\d+)\s+(\d+)/) && $flag_reverse==1){
	    $nb_mums_rev++;
	    
	    $tab_mums[$i]->{G0}->{sens}=1;
	    $tab_mums[$i]->{G0}->{init}=$1;
	    $tab_mums[$i]->{G0}->{long}=$3;
	    $tab_mums[$i]->{G0}->{fin}=$1+$3-1;
	    #modification due à l'option -c
	    $tab_mums[$i]->{G1}->{sens}=2;
	    $tab_mums[$i]->{G1}->{fin}=$2;
	    $tab_mums[$i]->{G1}->{long}=$3;
	    $tab_mums[$i]->{G1}->{init}=$2-$3+1;	  
	    $tab_mums[$i]->{ident}=$nb_mums_tot;  
	}
	$i++;
    }
    
}
close FIC;
    return \@tab_mums;
}


=head2 DoubleChev
  
  Description:
  
    Recherche les MUMs double chevauchants: elimination du MUM

   
  Arguments:
  
    $rtab_in  	Structure de donnees complexe contenant toutes les informations des MUMs restant apres
    			la suppression des inclusions dans les deux genomes.
    'G0' ou 'G1' selon le génome à traiter
  Returns:
  
    $rtab_out		Tableau contenant la reference vers la structure complexe modifie et la flag de savoir si, lors de la comparaison deux a deux des MUMs, au moins
    			un chevauchement a ete traite et donc de savoir si il va falloir relancer une autre
			fois la fonction.
    
=cut 
sub DoubleChev
{

	my $rtab_in=shift;
	my $gen=shift;
	my $cptr_del=0;
	my @tab_in=@$rtab_in;
	my $continue=0;
	my @baggins=();
	for (my $i=$#tab_in-1;$i>0;$i--){
	    if (($tab_in[$i-1]->{$gen}->{fin})>=($tab_in[$i+1]->{$gen}->{init}-1)&&($tab_in[$i]->{$gen}->{fin}<=$tab_in[$i+1]->{$gen}->{fin})){
		#cas du mum inclus dans ses deux voisins
		push (@baggins, $tab_in[$i]);
		$i--;
	    }
	}
	return \@baggins;
}

=head2 Remove_symetrically
  
  Description:
  
    A partir des deux listes de MUM inclus détectés, on merge en une seule, 
    La liste de mum est convertie en hash (cle=ident)
    On fait ensuite des deletions du hash, sur les cles correspondant 
    aux elements du merge

  Arguments:
  
    $rtab_in  	Structure de donnees complexe contenant toutes les informations des MUMs restant apres
    			la suppression des inclusions dans les deux genomes.
   $r_baggins, les sacs de mum à enlever, sur G0 et sur G1
  Returns:
  
    @tab_new		Tableau contenant la structure complexe modifie 

    
=cut 
sub Remove_symetrically
{
    my $rtab_in=shift;
    my $r_baggins_G0=shift;
    my $r_baggins_G1=shift;
    my @baggins_G0=@$r_baggins_G0;
    my @baggins_G1=@$r_baggins_G1;
    my @tab_in=@$rtab_in;
    my (%merge, %dic_in, @all, $l,$l1, $l2, $i, $j, $k, $key, $key_2_delete, @tab_new, $id);
    for ($i=0;$i <= $#baggins_G0;$i++){
	$id=$baggins_G0[$i]->{ident};
	$merge{$id}=$baggins_G0[$i];
    }
    @all=keys(%merge);
    my $l_interm=scalar(@all);
    for ($i=0;$i <= $#baggins_G1;$i++){
	$id=$baggins_G1[$i]->{ident};
	$merge{$id}=$baggins_G1[$i];
    }
    @all=keys(%merge);
    $l= scalar(@all);
    $l1= scalar(@baggins_G0);
    $l2=scalar(@baggins_G1);

    for ($k=0; $k<= $#tab_in; $k++){
	$key= $tab_in[$k]->{ident};
	$dic_in{$key}=$tab_in[$k];
    }
    for ($j=0; $j<= $#all; $j++){
	$key_2_delete= $all[$j];
	delete $dic_in{$key_2_delete};
    }
    @tab_new=values (%dic_in);
    return \@tab_new;
}

=head2 RechChev_trimbyend
  
  Description:
  
    Recherche les MUMs chevauchants sur le génome G0 et coupe la zone chevauchante par la fin du MUM sur G0, tout
    en repercutant cette diminution de taille sur le MUM correspondant dans l autre genome.
    La recherche se fait deux a deux en partant des MUMs de position initale sur le genome G0 la plus haute
    
   
  Arguments:
  
    $rtab_in  	Structure de donnees complexe contenant toutes les informations des MUMs restant apres
    			la suppression des inclusions dans les deux genomes.
    $gen        Génome à traiter vaut 'G0' ou 'G1'
			
  Returns:
  
    $rtab_new	Référence vers le tableau contenant la reference vers la structure complexe modifiée
    
=cut   
sub RechChev_trimbyend
{
	my $rtab_in=shift;
	my $gen=shift;
	my $cptr_trim=0;
	my @tab_in=@$rtab_in;
	my $tmp_len_ch=0;
	
	for (my $i=$#tab_in-1;$i>=0;$i--){

	    
	    if(($tab_in[$i]->{$gen}->{fin})>=($tab_in[$i+1]->{$gen}->{init})){
		#cas du mum chevauchant mais non inclus dans ses deux voisins	
		$cptr_trim++;
		
		if ($gen eq 'G0'){
		#trim sur G0
		    $tmp_len_ch=($tab_in[$i]->{G0}->{fin})-($tab_in[$i+1]->{G0}->{init})+1;
		    #on rogne par la fin
		    ($tab_in[$i]->{G0}->{fin})=($tab_in[$i+1]->{G0}->{init})-1;
		    $tab_in[$i]->{G0}->{long}=($tab_in[$i]->{G0}->{fin})-($tab_in[$i]->{G0}->{init})+1;

		    if($tab_in[$i]->{G1}->{sens}==1){
			$tab_in[$i]->{G1}->{fin}=($tab_in[$i]->{G1}->{fin})-$tmp_len_ch;
			$tab_in[$i]->{G1}->{long}=($tab_in[$i]->{G1}->{fin})-($tab_in[$i]->{G1}->{init})+1;

		    }else{

			$tab_in[$i]->{G1}->{init}=($tab_in[$i]->{G1}->{init})+$tmp_len_ch;
			$tab_in[$i]->{G1}->{long}=($tab_in[$i]->{G1}->{fin})-($tab_in[$i]->{G1}->{init})+1;
		    }										
		}

		else {
		# trim sur G1    			
		    $tmp_len_ch=($tab_in[$i]->{G1}->{fin})-($tab_in[$i+1]->{G1}->{init})+1;
		    ($tab_in[$i]->{G1}->{fin})=($tab_in[$i+1]->{G1}->{init})-1;
		    ($tab_in[$i]->{G1}->{long})=($tab_in[$i]->{G1}->{fin})-($tab_in[$i]->{G1}->{init})+1;

		    if($tab_in[$i]->{G1}->{sens}==1){
			($tab_in[$i]->{G0}->{fin})=($tab_in[$i]->{G0}->{fin})-$tmp_len_ch;
			($tab_in[$i]->{G0}->{long})=($tab_in[$i]->{G0}->{fin})-($tab_in[$i]->{G0}->{init})+1;

		    }else{
			($tab_in[$i]->{G0}->{init})=($tab_in[$i]->{G0}->{init})+$tmp_len_ch;
			($tab_in[$i]->{G0}->{long})=($tab_in[$i]->{G0}->{fin})-($tab_in[$i]->{G0}->{init})+1;
		    }		

		}
	    }
	}

	return $rtab_in;
}


=head2 RechInc
  
  Description:
  
    Recherche les MUMs inclus dans d autres MUMs sur le génome de choix et les supprime de la liste des MUMs
    La recherche se fait deux a deux en partant des MUMs de position initale sur le genome G0 la plus haute
    et de manière itérative tant qu une inclusion est detectee.
    
  Arguments:
  
    $rtab_in  	Structure de donnees complexe contenant toutes les informations des MUMs contenues
    			dans la sortie de Mummer.
    'G0' ou 'G1' selon le genome a traiter			
  Returns:
  
    $rtab_out		Tableau contenant la reference vers la nouvelle structure complexe et la flag permettant 
                        de savoir si, lors de la comparaison deux a deux des MUMs, au moins
    			une inclusion a ete traitee et donc de savoir si il va falloir relancer une autre
			fois la fonction.
    
=cut 

sub RechInc
{
	my $rtab_in=shift;
	my $gen=shift;
	my @tab_in=@$rtab_in;
	my @collection=();
	my $i=0;
	my $j;
	while ($i<$#tab_in)
	{
	    if( (($tab_in[$i]->{$gen}->{fin}) >= ($tab_in[$i+1]->{$gen}->{fin})))
	    {
		push(@collection, $tab_in[$i+1]);
		$j=$i;
		$j++;
		while (($j < $#tab_in) && (($tab_in[$i]->{$gen}->{fin}) >= ($tab_in[$j+1]->{$gen}->{fin})) )
		{
		    push(@collection, $tab_in[$j+1]);
		    $j++;
		}
		$i=$j+1;
	    }
	    else 
	    {
		$i++;
	    }
	}
		    
	return \@collection;
}
	    

sub tell_me
{
my $rtab_in=shift;
my @tab_in=@$rtab_in;
my $somme_mums_nn_ch=0;
for (my $i=0;$i<scalar(@tab_in);$i++){
     $somme_mums_nn_ch+=$tab_in[$i]->{G0}->{long};
 }
my $nb_mums_nn_ch=scalar(@tab_in);
my $distance;
$distance=1-($somme_mums_nn_ch/(($gen1_len+$gen2_len)/2));

return $distance;
}

sub order_G0
{
$a->{G0}->{init} <=> $b->{G0}->{init}
		 ||
$b->{G0}->{long} <=> $a->{G0}->{long}
                 ||
$a->{G1}->{init} <=> $b->{G1}->{init}

}

sub order_G1
{
$a->{G1}->{init} <=> $b->{G1}->{init}
		 ||
$b->{G1}->{long} <=> $a->{G1}->{long}
		 ||
$a->{G0}->{init} <=> $b->{G0}->{init}
}


sub pretreatment
{
    my $rtab_in=shift;
    my @tab_in_trie_G1;my @tab_in_trie_G0;
    my $collec_G0; my $collec_G1;	
    my $rtab_new;
	
    #Pre-traitement des inclusions sur G0 et G1
    @tab_in_trie_G0=sort order_G0 @$rtab_in;
    $collec_G0=RechInc(\@tab_in_trie_G0, 'G0');
    @tab_in_trie_G1=sort order_G1 @tab_in_trie_G0;
    $collec_G1=RechInc(\@tab_in_trie_G1, 'G1');
    $rtab_new=Remove_symetrically($rtab_in, $collec_G0, $collec_G1);
    return $rtab_new;
}

sub treat_dc
{
    my $rtab_in=shift;
    my @tab_noincl_og0; my @tab_noincl_og1;
    my $baggins_G0; my $baggins_G1;
    my $rtab_new;

    #Elimination doubles chevauchants sur G0 et G1
    @tab_noincl_og0=sort order_G0 @$rtab_in;
    $baggins_G0=DoubleChev(\@tab_noincl_og0, 'G0');
    @tab_noincl_og1=sort order_G1 @$rtab_in;
    $baggins_G1=DoubleChev(\@tab_noincl_og1, 'G1');
    $rtab_new=Remove_symetrically($rtab_in, $baggins_G0, $baggins_G1);
				
    return $rtab_new;
}

sub treat_chG0
{
    my $rtab_in=shift;
    my @tab_nodc_og0;
    my $rtab_new;

    #Traite chevauchants simples de G0
    @tab_nodc_og0=sort order_G0 @$rtab_in;
    $rtab_new=RechChev_trimbyend(\@tab_nodc_og0, 'G0');
				
    return $rtab_new;
}

sub treat_chG0_d
{
    my $rtab_in=shift;
    my @tab_nodc_og0;
    my $rtab_new;

    #Traite chevauchants simples de G0
    @tab_nodc_og0=sort order_G0 @$rtab_in;
    $rtab_new=RechChev_trimbyend(\@tab_nodc_og0, 'G0');
				
    my $dist=tell_me($rtab_new);
    return $dist;
}

sub preG1_postG0
{
    my $rtab_in=shift;
    my @tab_in_trie_G1;
    my @collec_G0 =(); my $collec_G1; 
    my $rtab_new;
		
    #Pre-traitement G1, post G0
    @tab_in_trie_G1=sort order_G1 @$rtab_in;
    $collec_G1=RechInc(\@tab_in_trie_G1, 'G1');
    $rtab_new=Remove_symetrically($rtab_in, \@collec_G0, $collec_G1);
		
    return $rtab_new;
}

sub preG0_postG1
{
    my $rtab_in=shift;
    my @tab_in_trie_G0;
    my @collec_G1 =(); my $collec_G0; 
    my $rtab_new;
		
    #Pre-traitement G1, post G0
    @tab_in_trie_G0=sort order_G0 @$rtab_in;
    $collec_G0=RechInc(\@tab_in_trie_G0, 'G0');
    $rtab_new=Remove_symetrically($rtab_in, \@collec_G1, $collec_G0);
		
    return $rtab_new;
}

sub treat_chG1
{
    my $rtab_in=shift;
    my @tab_nodc_og1;
    my $rtab_new;
	
    #Traite chevauchants simples de G1
    @tab_nodc_og1=sort order_G1 @$rtab_in;
    $rtab_new=RechChev_trimbyend(\@tab_nodc_og1,'G1');

    my $dist=tell_me($rtab_new);
    return $dist;
}

sub treat_chG1_nod
{
    my $rtab_in=shift;
    my @tab_nodc_og1;
    my $rtab_new;
	
    #Traite chevauchants simples de G1
    @tab_nodc_og1=sort order_G1 @$rtab_in;
    $rtab_new=RechChev_trimbyend(\@tab_nodc_og1,'G1');

    return $rtab_new;
}

sub Give_MUMI1
{
    my $file_mums=shift;
    my $out=Read_mums($file_mums);
    $out= pretreatment($out);
    $out=treat_dc($out);
    $out=treat_chG0($out);
    $out=preG1_postG0($out);
    $out=treat_chG1($out);

    return $out;
}

sub Give_MUMI2
{
    my $file_mums=shift;
    my $out=Read_mums($file_mums);
    $out= pretreatment($out);
    $out=treat_dc($out);
    $out=treat_chG1_nod($out);
    $out=preG0_postG1($out);
    $out=treat_chG0_d($out);
    return $out;
}



my $mumi1=Give_MUMI1($ficName);
my $mumi2=Give_MUMI2($ficName);
my $AV_MUMI=($mumi1+$mumi2)/2;
print "$ficName $AV_MUMI\n";

