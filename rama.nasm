section .data
    alfa_phi dq -57.8
    alfa_psi dq -47.0
    beta_phi dq -119.0
    beta_psi dq 113.0
    half dq 0.5

section .text
    global rama
    rama:
        push rbp
        mov rbp, rsp

        vbroadcastsd ymm8, qword[alfa_phi]
        vbroadcastsd ymm9, qword[alfa_psi]
        vbroadcastsd ymm10, qword[beta_phi]
        vbroadcastsd ymm11, qword[beta_psi]
        vbroadcastsd ymm12, qword[half]

        vxorpd ymm0, ymm0, ymm0
        xor r8, r8

.loop_start:
        cmp r8, rdi
        jge .loop_end

        vmovupd ymm1, [rsi+r8*8]
        vmovupd ymm2, [rdx+r8*8]

        ; Calcolo distanza alfa
        vsubpd ymm3, ymm1, ymm8
        vmulpd ymm3, ymm3, ymm3

        vsubpd ymm4, ymm2, ymm9
        vmulpd ymm4, ymm4, ymm4

        vaddpd ymm3, ymm3, ymm4
        vsqrtpd ymm3, ymm3

        ; Calcolo distanza beta
        vsubpd ymm5, ymm1, ymm10
        vmulpd ymm5, ymm5, ymm5

        vsubpd ymm6, ymm2, ymm11
        vmulpd ymm6, ymm6, ymm6

        vaddpd ymm5, ymm5, ymm6
        vsqrtpd ymm5, ymm5

        ; Calcolo energia
        vminpd ymm7, ymm3, ymm5
        vmulpd ymm7, ymm7, ymm12
        vaddpd ymm0, ymm0, ymm7

        add r8, 4
        jmp .loop_start

.loop_end:
        vextractf128 xmm1, ymm0, 1
        vaddpd xmm0, xmm0, xmm1
        vhaddpd xmm0, xmm0, xmm0

        mov rsp, rbp
        pop rbp
        ret

