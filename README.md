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
