Algoritmi su implementirani u programskom jeziku C i prevođenje se vrši naredbom:

gcc ime_datoteke.c -o ime_datoteke
npr. gcc risotto.c -o risotto

Za izvršavanje potrebno je proslediti argumente komandne liste.
Jedini argument koji je potreban svim algoritmima jeste datoteka sa ulaznim sekvencama.
Ovo je obična .txt datoteka koja ima određen broj redova DNK sekvenci.
Pored njega postoji jedan argument koji je opcioni, a to je azbuka.
Ukoliko se ne navede, podrazumevana je A, C, G i T, ali se to može promeniti prosleđivanjem .txt fajla sa drugim karakterima.

Ostali argumenti se razlikuju tako da će biti opisani za svaki algoritam.


Potrebni argumenti za algoritam grube sile su: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]
gcc brute_force.c -o brute_force
./brute_force 13 3 ulazne_sekvence.txt azbuka.txt


Potrebni argumenti za algoritam genomski su: <duzina_motiva> <datoteka_sa_sekvencama> <verovatnoca_mutacija [0,1]> <broj_iteracija> [<datoteka_sa_azbukom>]
gcc ga.c -o ga
./ga 13 ulazne_sekvence.txt 0.2 100 azbuka.txt


Potrebni argumenti za algoritam MITRA su: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]
gcc mitra.c -o mitra
./mitra 13 3 ulazne_sekvence.txt azbuka.txt

Potrebni argumenti za algoritam PMS5 su: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> <datoteka_ilr_tabela> [<datoteka_sa_azbukom>]
gcc ga.c -o ga
./ga 13 3 ulazne_sekvence.txt ilr_tabela.txt azbuka.txt

ilr_tabela je .txt datoteka koja u sebi sadrži rešenje linearnih kombinacija za određene vrednosti dužine motiva i broja mutacija. Sadrži osmodimenzioni niz i bool vrednosti npr. 
0 0 1 1 0 0 2 2 true.


Potrebni argumenti za algoritam slučajnu projekciju i EM su: <duzina_motiva> <broj_proj> <s-filtriranje_korpi> <iteracije> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]
gcc random_projection_and_em.c -o random_projection_and_em
./random_projection_and_em 13 6 4 50 ulazne_sekvence.txt azbuka.txt


Potrebni argumenti za algoritam Risotto su: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]
gcc risotto.c -o risotto
./risotto 13 3 ulazne_sekvence.txt azbuka.txt


Potrebni argumenti za glasački algoritam su: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]
gcc voting.c -o voting
./voting 13 3 ulazne_sekvence.txt azbuka.txt

Potrebni argumenti za algoritam Winnower su: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> <k-uslov odsecanja> [<datoteka_sa_azbukom>]
gcc winnower.c -o winnower
./winnower 13 3 ulazne_sekvence.txt 3 azbuka.txt

