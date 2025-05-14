# Progetto Architetture Avanzate

Questo progetto è stato sviluppato come parte del corso di Architetture Avanzate, con l'obiettivo di implementare un algoritmo di ottimizzazione basato sul Simulated Annealing per la modellazione di strutture proteiche.  
Il focus principale è sull'efficienza computazionale, sfruttando l'assembly per ottimizzare le prestazioni delle operazioni critiche.

## Struttura del Progetto

- **`pst32c.c` e `pst64c.c`**: Implementazioni in C per architetture a 32 e 64 bit, rispettivamente.
- **`rama.nasm`**: Routine in assembly NASM per l'ottimizzazione di calcoli specifici legati all'energia conformazionale.
- **`README.md`**: Documentazione del progetto e note di sviluppo.

## Funzionalità Implementate

- **Funzione di rotazione**: Implementata e testata; funzionante ma con margini di ottimizzazione in assembly.
- **Prodotto matriciale**: Implementato nel modulo "backbone"; è necessario verificare la correttezza del passaggio dei parametri.
- **Calcolo delle energie**:
  - `packing_energy`: Restituisce un valore molto alto; da verificare.
  - Altre energie (escluse `rama-energy`): Restituiscono 0; è necessario approfondire.
  - `energia_tot`: Risulta in un numero molto grande; richiede ulteriori analisi.
- **Algoritmo di Simulated Annealing**: Implementato; durante i test con input casuali si verifica un segmentation fault, probabilmente legato alla gestione del vettore `coords`.

## Note di Sviluppo

- Il file di input deve essere un binario codificato in un formato specifico.
- Sono stati identificati problemi nel vettore `coords`, con valori anomali che portano a segmentation fault.
- È necessario verificare la gestione della memoria e la correttezza dei dati durante l'elaborazione.

## Requisiti

- Compilatore C compatibile con lo standard C99 o superiore.
- Assembler NASM per la compilazione del file `rama.nasm`.
- Sistema operativo Linux o compatibile.

## Compilazione

Per compilare il progetto:

```bash
nasm -f elf64 rama.nasm -o rama.o
gcc -m64 pst64c.c rama.o -o progetto
```

Assicurarsi di utilizzare le opzioni corrette per l'architettura target (32 o 64 bit).

## Esecuzione

Dopo la compilazione, eseguire il programma passando il file binario di input:

```bash
./progetto input.bin
```

## Contributi

Il progetto è aperto a contributi. Si prega di aprire una issue per discutere modifiche o miglioramenti.

## Autori

- [Andrea Farfaglia](https://github.com/andreaahahah)
- [Lorenzo Manna](https://github.com/lmann97)
- [Angelica Porco](https://github.com/AngeP02)
