#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

long double factorial(unsigned n){
	long double f=1;
	while(n>0){f*=n--;}
	return f;
}

long double combination(unsigned k,unsigned n){
	long double f=(factorial(n) / (factorial(k) * factorial(n-k)));
	return f;
}


long double combination_eff(unsigned k,unsigned n){
	long double num=1;
	if (k<(n/2)){k=n-k;}
	int n_2=n;
	while (n_2>k){num*=n_2--;}
	long double f= num / factorial(n-k);
	return f;
}

long double proba_n(unsigned n,unsigned k, long double proba){
	long double result=combination_eff(k,n) * powl(proba,k) * powl((1-proba),(n-k)); // New way more efficient to compute combination
	return result;
}


long double proba_more_than(int n,int k, long double proba){
	long double result=0.0;
	while(k<=n) {
		result+=combination_eff(k,n) * powl(proba,k) * powl((1-proba),(n-k));
		k++;
	}
	return result;
}


long double proba_less_than(int n,int k, long double proba){
	long double result=0.0;
	while(k>=0) {
		result+=combination_eff(k,n) * powl(proba,k) * powl((1-proba),(n-k));
		k--;
	}
	return result;
}

int get_th(int size_window,long double threshold, long double proba){
	int th_nb_gene=size_window+1;
	long double p_t=0.0;
// 	printf("starting at %d / with proba %LE\n",th_nb_gene,proba);
	while(p_t<=threshold && th_nb_gene>0){
		th_nb_gene--;
		p_t = p_t + proba_n(size_window,th_nb_gene,proba);
// 		printf("\tp(x>=%d) = %LE\n",th_nb_gene,p_t);
	}
	return th_nb_gene;
}


int get_th_less(int size_window,long double threshold, long double proba){
	int th_nb_gene=-1;
	long double p_t=0.0;
	while(p_t<=threshold && th_nb_gene<size_window){
		th_nb_gene++;
		p_t = p_t + proba_n(size_window,th_nb_gene,proba);
	}
	return th_nb_gene;
}


int is_local_maximum(int start,int size,int type, int p_nb_genes, int p_max,double ***store){
	int i,j,result=1,hood=5; // What's the size of the 'hood 'bro ?
// 	printf("we'll look for 5 around %d and for 5 around %d with a maximum index of %d and %d\n",start,size,p_nb_genes,p_max);
	for (i=start-hood;i<=start+hood;i++){
		for(j=size-hood;j<=size+hood;j++){
// 			printf("-- Looking at %d %d\n",i,j);
			if (i>=0 && j>=0 && i<p_nb_genes && j<=p_max){ // should be enough to get if there is a value in the table
// 				printf("-- Really Looking at %d %d\n",i,j);
				if (store[i][j][type]>store[start][size][type]){
					result=0;
// 					i=start+hood+1;j=size+hood+1;
				}
			}
		}
	}
	return result;
}


long double log10perso(long double x){
	return log(x)/log(10);
}


int main(int argc, char *argv[])
{
// 	printf( "I am alive!  Beware.\n" );
	FILE *ifp, *reffile;
	char* refFilename=argv[1];char* inputFilename=argv[2];char* outputFilename=argv[3];
	reffile=fopen(refFilename,"r");
	int nb_genes=0,phage=0,pfam=0,unch=0,size=0,strand=0,hallmark=0,i=0,noncaudo=0;
	float f_size=0.0;
	long double p_phage=0.0,p_pfam=0.0,p_unch=0.0,p_strand=0.0,p_noncaudo=0.0;
	if (reffile == NULL) {
		fprintf(stderr, "Can't open input file %s\n",refFilename);
		exit(1);
	}
	while (fscanf(reffile,"%Lf %Lf %Lf %Lf %f %Lf", &p_phage, &p_pfam, &p_unch, &p_strand, &f_size, &p_noncaudo) == 6) {}
	printf("refs => %LE %LE %LE %LE %f %LE\n", p_phage, p_pfam, p_unch, p_strand, f_size, p_noncaudo);
	fclose(reffile);
	ifp = fopen(inputFilename, "r");
	if (ifp == NULL) {
		fprintf(stderr, "Can't open input file %s!\n",inputFilename);
		exit(1);
	}
	if (fscanf(ifp, "%d", &nb_genes) == 1){
// 		printf("%d genes\n",nb_genes);
	}
	// Alloc memory for gene tables
	int t_phage[nb_genes],t_pfam[nb_genes],t_unch[nb_genes], t_size[nb_genes],t_strand[nb_genes],t_hallmark[nb_genes],t_noncaudo[nb_genes];
	while (fscanf(ifp,"%d %d %d %d %d %d %d", &phage, &noncaudo, &pfam, &unch, &size, &strand, &hallmark) == 7) {
// 		printf("gene %d => %d %d %d %d %d %d %d\n", i, phage, noncaudo, pfam, unch, size, strand, hallmark);
		t_phage[i]=phage;
		t_noncaudo[i]=noncaudo;
		t_pfam[i]=pfam;
		t_unch[i]=unch;
		t_size[i]=size;
		t_strand[i]=strand;
		t_hallmark[i]=hallmark;
		i++;
	}
	fclose(ifp);
	if (nb_genes!=i){
		printf("Houston we got a problem !!!!!! : we had %d genes and we count %d lines\n",nb_genes,i);
		exit(1);
	}
// 	// set up sliding windows
	int min=10,max=100;
	if (min>nb_genes){min=nb_genes;}
	if (max>nb_genes){max=nb_genes;}
// 	// how many sliding windows will we have ?
	int k=0,j=0,max_g=0,c_phage=0,c_pfam=0,pred_nb_s_w=0,t=0,th_nb_gene=0;
	for (k=min;k<=max;k++){
		pred_nb_s_w+=nb_genes-k+1;
	}
// 	printf("Predicting %d sliding windows\n",pred_nb_s_w);
	// computing the threshold for each size of sliding window
// 	printf("Trying to allocate the memory 1\n");
	long double th=0.01/pred_nb_s_w,p_t=0.0;
	// alloc memory for score matrix for the 6 metrics
	double ***store=malloc(nb_genes*sizeof(double **));
	if (store==NULL){printf("out of memory\n");exit(1);}
	for(i=0; i < nb_genes; i++){
		store[i] = malloc(max * sizeof(double *));
		if(store[i] == NULL){printf("out of memory\n");exit(1);}
		for (j=0;j<=max;j++){
			store[i][j] = malloc(6 * sizeof(double ));
			if(store[i][j] == NULL){printf("out of memory\n");exit(1);}
			for (k=0;k<6;k++){store[i][j][k]=0;}
		}
	}
// 	printf("Memory Allocated and Initialized for %d %d 5\n",nb_genes,max);
	int store_h[nb_genes][max];
	int n_phage=0,n_pfam=0,n_short=0,n_switch=0,n_unch=0,n_hallmark=0,n_noncaudo=0;
	printf("For this contig we'll have %d sliding windows (= nb of comparison)\n",pred_nb_s_w);
	for (k=max;k>=min;k--){
		int th_phage=k,th_pfam=k,th_size=k,th_unch=k,th_strand=k,th_noncaudo=k;
		// we get all thresholds
		th_phage=get_th(k,th,p_phage);
// 		printf("For window size %d, you will need at least %d phage genes to be significant\n",k,th_phage);
		th_pfam=get_th_less(k,th,p_pfam);
		th_unch=get_th(k,th,p_unch);
// 		printf("For window size %d, you will need at least %d uncharacterized genes to be significant\n",k,th_unch);
		th_size=get_th(k,th,0.1);
		th_strand=get_th_less(k,th,p_strand);
		th_noncaudo=get_th(k,th,p_noncaudo);
// 		printf("For window size %d, you will need at least %d noncaudo genes to be significant\n",k,th_noncaudo);
// 		printf("////// Sliding window of %d genes -> th %d\n",k,th_phage);
		// For all the sliding windows of this size, we count and compute and store the significativity value if > sig
		for (i=0;i<(nb_genes-k+1);i++){
			n_phage=0;n_pfam=0;n_unch=0;n_short=0;n_switch=0;n_hallmark=0;n_noncaudo=0;
// 			// Counting 
			for (j=i;j<(i+k);j++){
				n_phage+=t_phage[j];
// 				printf("Adding %d to the number of phage genes (%d)\n",t_phage[j],j);
				n_pfam+=t_pfam[j];
				n_unch+=t_unch[j];
				n_short+=t_size[j];
				n_switch+=t_strand[j];
				n_hallmark+=t_hallmark[j];
				n_noncaudo+=t_noncaudo[j];
			}
			unsigned tag=0;
// 			// If above thresholds
			if (n_phage>th_phage){
// 				// Calculate and store significativity
				store[i][k][0]=-1*log10(proba_more_than(k,n_phage,p_phage)*pred_nb_s_w);tag=1;
// 				printf("Phage => %d is beyond the threshold %d, so we compute its significativity %E, that we store in %d, %d, 0\n",n_phage,th_phage,store[i][k][0],i,k);
			}
			if (n_pfam<th_pfam){
				// Calculate and store significativity
				store[i][k][1]=-1*log10(proba_less_than(k,n_pfam,p_pfam)*pred_nb_s_w);tag=1;
// 				printf("Pfam => %d is below the threshold %d, so we compute its significativity %E, that we store in %d, %d, 1\n",n_pfam,th_pfam,store[i][k][1],i,k);
			}
			if (n_unch>th_unch){
// 				// Calculate and store significativity
				store[i][k][2]=-1*log10(proba_more_than(k,n_unch,p_unch)*pred_nb_s_w);tag=1;
// 				printf("Unch => %d is beyond the threshold %d, so we compute its significativity %E, that we store in %d, %d, 2\n",n_unch,th_unch,store[i][k][2],i,k);
			}
			if (n_short>th_size){
// 				// Calculate and store significativity
				store[i][k][3]=-1*log10(proba_more_than(k,n_short,0.1)*pred_nb_s_w);tag=1;
// 				printf("Short => %d is beyond the threshold %d, so we compute its significativity %E, that we store in %d, %d, 3\n",n_short,th_size,store[i][k][3],i,k);
			}
			if (n_switch<th_strand){
				// Calculate and store significativity
				store[i][k][4]=-1*log10(proba_less_than(k,n_switch,p_strand)*pred_nb_s_w);tag=1;
// 				printf("Switch => %d is below the threshold %d, so we compute its significativity %E, that we store in %d, %d, 4\n",n_switch,th_strand,store[i][k][4],i,k);
			}
			if (n_noncaudo>th_noncaudo){
// 				// Calculate and store significativity
				store[i][k][5]=-1*log10(proba_more_than(k,n_noncaudo,p_noncaudo)*pred_nb_s_w);tag=1;
// 				printf("Phage => %d is beyond the threshold %d, so we compute its significativity %E, that we store in %d, %d, 0\n",n_phage,th_phage,store[i][k][0],i,k);
			}
			if (tag==1){store_h[i][k]=n_hallmark;}
		}
	}
	// We look for local maxima and export the results
	FILE *ofp;
	ofp = fopen(outputFilename, "w");
	if (ofp == NULL) {
		fprintf(stderr, "Can't open output file %s!\n",outputFilename);
		exit(1);
	}
	for (k=max;k>=min;k--){
		for (i=0;i<(nb_genes-k+1);i++){
			for (j=0;j<6;j++){
				if (store[i][k][j] != 0.0){ // the stored value is not null
// 					printf("potential local maximum %d %d %d %E %d\n",i,k,j,store[i][k][j],store_h[i][k]);
					if (is_local_maximum(i,k,j,nb_genes-1,max,store)==1){ // and is a local maxima
						// so we print it, with the nb_hallmark (start / window size / type / sig / nb_hallmark)
// 						printf("local maximum ! %d %d %d %E %d\n",i,k,j,store[i][k][j],store_h[i][k]);
						// i - start gene / k - sliding window size / j - proof typ (0 - phage / 1 - pfam / 2 - unch / 3 - size / 4 - strand)
						fprintf(ofp, "%d\t%d\t%d\t%.14lf\t%d\n",i,k,j,store[i][k][j],store_h[i][k]);
					}
				}
			}
		}
	}
	fclose(ofp);
	printf("done");
	// We export the results
	return 0;
}
