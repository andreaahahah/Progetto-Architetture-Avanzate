implementata funzione rotation, testata e funzionante ma da ottimizzare in assembly

TODO
implementare prodotto matriciale in backbone
verificare la correttezza della funzione con il passaggio dei parametri


19/12
implementata backbone ed inoltre implementate le varie energie,ci danno dei valori ma va verificato se sono corretti(packing_energy ha un valore molto alto,le altre due(non considerando rama-energy) danno 0 come valore->risultato:energia_tot un numero molto grande.

21/12
implementato l'algoritmo di simulated annealing,provato inoltre a fare una prima fase di testing con input casuali ma va in segmentation fault,da rivedere.


27/12
risolto il problema del file da passare nella sequenza
deve essere un binario codificato in una certa maniera
facendo delle print abbiamo individuato il problema nel vettore coords
facendo delle print si vede come i valori di coords siano sballati e poi va verificato perchè anche cambiando la dimensione non va in seg.fault
facendo delle print risultano valori sballati anche in newv,però lui prende solo la funzione rotation e le distanza ca_n ecc. che sono già definite a priori;quindi va verificata anche rotation perchè potrebbe esserre l'errore a monte che ha scatenato le  coordinate errate.


30/12

i valori degli angoli phi e psi danno valori molto elevati nei file di output
abbiamo risolto i problemi dei valori sballati in rotation e backbone usando al posto delle approssimazioni lineari 
basate su serie di taylor del sen e cos fatte da noi le funzioni sin e cos della libreria math.h
il problema ora sembra essre nelle energie perchè elec ad esmpio rimane sempre a 0
e packing un valore molto alto 
abbiamo ristretto il problema al fatto di come prendevamo i valori ascii relativi alla posizione degli amminoacidi nella sequenza,quindi abbiamo provato a ristampare e alcuni valori comunque non rientrano nel range 
ora le prime tre energie comprese elec sembrano essere corrette, ma packing stampa sempre valori molto alti ;
abbiamo controllato il simulated annealing e sembra che arriviamo ad un minimo globale quindi sembra corretta


03/01



Risolto il probelma degli angoli facendo la normalizzazione;ora il range degli angoli sembra essere corretto.



10/01
  *********************************************************************************************************************************************************************
effettuato ricevimento dal prof;
primo problema risolto;doppia generaziome casuale degli angoli phi e psi fatti già da prima da gen_rnd_mat che andava a sovrapporre i vecchi valori ed aggiornarli con altri senza trovarsi
secondo problema risolto:nella funzione rotation è stata tolta sqrt perchè è solo un prodotto scalare senza normalizzazione


terzo problema da aggiustare:
le norme in backbone di v1 v2 e v3 vanno fatte fuori dal for e non dentro il for;

**********************************************************************************************************************************************************************
correzioni da fare

aggiustare norma che va messa fuori dal for in backbone sia in v1 v2 e v3
togliere le print e le exit 0 rimaste
confrontare ris phi e psi
confrontare energie del prof e nostre
valori energie dei prof------->

rama_e: 9538.258065, hydro_e: 98.171527, elec_e: -0.540941, pack_e: 36599998.039925
ENERGIA INIZIALE: 10989586.647618
-------------------------------------->



20/01

risolto problema file pst ; ora i risultati corretti;
solo il file c ci dava una media di 0.6 sec
file a 64 bit
ottimizzata in assembly funzione rama_energy con code vectorization

con assembly:
angelica ->virtual machine ->circa 0.9 sec

lorenzo->circa 1.3 sec WSL(mini vm)

andrea->ambiente linux nativo con 8 core->circa 0.62 s 
//todo
fare qualche altra ottimizzazione

fare ottimizzazione rama in sse a 32

cont.20 gennaio

incominciata ottimizzazione pack
diversi errori risolti sintattici
loop infinito->vedere possibili errori

