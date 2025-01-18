%include "sseutils64.nasm"

section .data
    align 16
    alfa_phi dq -57.8
    alfa_psi dq -47.0
    beta_phi dq -119.0
    beta_psi dq 113.0
    half dq 0.5

section .bss
    res     resq 1      ; variabile per il risultato finale


section .text
global rama

rama:
   start
  ; sub rsp, 32                  ; riserva spazio per eventuali variabili locali
   mov ecx, edi                 ; interno n
   ;mov rax, rsi                 ; puntatore a phi
   ;vmovsd xmm0, qword[rax]      ; xmmo contiene puntatore a phi
   ;mov rax, rdx                 ; puntatore a psi   
   ;vmovsd xmm1,qword[rax]       ; xmm1 contiene puntatore a psi

   vxorpd ymm10, ymm10, ymm10 ; inizializza ymm10 (energia accumulata a zero)
 ; ------ fino a qua giusto ------

  
loop_unrolling:
        cmp ecx, 4              ; verifica se ci sono almeno 4 elementi da processare
        jl loop_end             ; se meno di 4, termina il loop
      
        ; carica 4 valori di phi e psi
        mov rax, rsi   
        vmovsd xmm0, [rax] 
      
        mov rax, rdx   
        vmovsd xmm1, [rax] 
        ; Calcolo alfa_dist e beta_dist
        vsubpd ymm2, ymm0, [alfa_phi]  ; |phi[i] - alfa_phi|
        vmulpd ymm2, ymm2, ymm2  ; (phi[i] - alfa_phi)^2     
        
        vsubpd ymm3, ymm1, [alfa_psi]  ; |psi[i] - alfa_psi|    
        vmulpd ymm3, ymm3, ymm3  ; (psi[i] - alfa_psi)^2
        

        vaddpd ymm2, ymm2, ymm3  ; somma dei quadrati
        vsqrtpd ymm2, ymm2       ; Calcola radice quadrata
        ; ymm2 contiene la radice quadrata della somma dei quadrati di alpha
        

        vsubpd ymm4, ymm0, [beta_phi]  ; |phi[i] - beta_phi|
        vmulpd ymm4, ymm4, ymm4  ; (phi[i] - beta_phi)^2  
    
        vsubpd ymm5, ymm1, [beta_psi]  ; |psi[i] - beta_psi|
        vmulpd ymm5, ymm5, ymm5  ; (psi[i] - beta_psi)^2

        vaddpd ymm4, ymm4, ymm5  ; somma dei quadrati
        vsqrtpd ymm4, ymm4       ; Calcola radice quadrata
        ; ymm4 contiene la radice quadrata della somma dei quadrati di beta

        
        vminpd ymm6, ymm2, ymm4  ; Trova il minimo di alfa e beta distanze


        
        vmulpd ymm6, ymm6, [half] ; moltiplica il minimo per 0.5
        
        vaddpd ymm10, ymm10, ymm6 ; accumula in ymm10        
        
        ;avanza i puntatori e decrementa la variabile i



        add rsi, 32             ; avanza phi di 4 elementi
        add rdx, 32             ; avanza psi di 4 elementi
        sub ecx, 4              ; decrementa la lunghezza residua
            
        jmp loop_unrolling         ; ripeti il ciclo

loop_end:
        ; somma orizzontale dei risultati accumulati
       
        vhaddpd ymm10, ymm10, ymm10     ; somma i valori in ymm10
        vhaddpd ymm10, ymm10, ymm10



        ;movsd [xmm2], ymm10
        vextractf128 xmm2, ymm10, 0
        vmovsd xmm0,xmm2
        ;movq xmm0, xmm2                  ; copia il risultato in r10
    
        ;salva in risultato in memoria
       ; movsd [res], xmm0


; ---- da qui in poi giusto ----

  ; add rsp, 32
   stop
  
