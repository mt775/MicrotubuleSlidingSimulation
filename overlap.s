	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 10, 15	sdk_version 10, 15, 6
	.globl	_overlap                ## -- Begin function overlap
	.p2align	4, 0x90
_overlap:                               ## @overlap
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$80, %rsp
	movl	%edi, -12(%rbp)
	movl	%esi, -16(%rbp)
	movl	-12(%rbp), %eax
	cmpl	-16(%rbp), %eax
	jne	LBB0_2
## %bb.1:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -8(%rbp)
	jmp	LBB0_11
LBB0_2:
	movl	-12(%rbp), %eax
	shll	$1, %eax
	subl	$1, %eax
	movl	%eax, -60(%rbp)
	movl	-12(%rbp), %eax
	shll	$1, %eax
	movl	%eax, -64(%rbp)
	movl	-16(%rbp), %eax
	shll	$1, %eax
	subl	$1, %eax
	movl	%eax, -68(%rbp)
	movl	-16(%rbp), %eax
	shll	$1, %eax
	movl	%eax, -72(%rbp)
	movl	-60(%rbp), %edi
	movl	-68(%rbp), %esi
	callq	_soverlap
	movsd	%xmm0, -32(%rbp)
	movl	-60(%rbp), %edi
	movl	-72(%rbp), %esi
	callq	_soverlap
	movsd	%xmm0, -40(%rbp)
	movl	-64(%rbp), %edi
	movl	-68(%rbp), %esi
	callq	_soverlap
	movsd	%xmm0, -48(%rbp)
	movl	-64(%rbp), %edi
	movl	-72(%rbp), %esi
	callq	_soverlap
	movq	_ovlp@GOTPCREL(%rip), %rcx
	movsd	%xmm0, -56(%rbp)
	movq	(%rcx), %rcx
	movslq	-60(%rbp), %rdx
	movq	(%rcx,%rdx,8), %rcx
	movslq	-68(%rbp), %rdx
	cmpl	$0, (%rcx,%rdx,4)
	jne	LBB0_4
## %bb.3:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -32(%rbp)
LBB0_4:
	movq	_ovlp@GOTPCREL(%rip), %rax
	movq	(%rax), %rax
	movslq	-60(%rbp), %rcx
	movq	(%rax,%rcx,8), %rax
	movslq	-72(%rbp), %rcx
	cmpl	$0, (%rax,%rcx,4)
	jne	LBB0_6
## %bb.5:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -40(%rbp)
LBB0_6:
	movq	_ovlp@GOTPCREL(%rip), %rax
	movq	(%rax), %rax
	movslq	-64(%rbp), %rcx
	movq	(%rax,%rcx,8), %rax
	movslq	-68(%rbp), %rcx
	cmpl	$0, (%rax,%rcx,4)
	jne	LBB0_8
## %bb.7:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -48(%rbp)
LBB0_8:
	movq	_ovlp@GOTPCREL(%rip), %rax
	movq	(%rax), %rax
	movslq	-64(%rbp), %rcx
	movq	(%rax,%rcx,8), %rax
	movslq	-72(%rbp), %rcx
	cmpl	$0, (%rax,%rcx,4)
	jne	LBB0_10
## %bb.9:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -56(%rbp)
LBB0_10:
	movsd	-32(%rbp), %xmm0        ## xmm0 = mem[0],zero
	addsd	-40(%rbp), %xmm0
	addsd	-48(%rbp), %xmm0
	addsd	-56(%rbp), %xmm0
	movsd	%xmm0, -24(%rbp)
	movsd	-24(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movsd	%xmm0, -8(%rbp)
LBB0_11:
	movsd	-8(%rbp), %xmm0         ## xmm0 = mem[0],zero
	addq	$80, %rsp
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.section	__TEXT,__literal8,8byte_literals
	.p2align	3               ## -- Begin function soverlap
LCPI1_0:
	.quad	4602678819172646912     ## double 0.5
LCPI1_2:
	.quad	4562254508917369340     ## double 0.001
	.section	__TEXT,__literal4,4byte_literals
	.p2align	2
LCPI1_1:
	.long	1073741824              ## float 2
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_soverlap
	.p2align	4, 0x90
_soverlap:                              ## @soverlap
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$96, %rsp
	movl	%edi, -12(%rbp)
	movl	%esi, -16(%rbp)
	movl	-12(%rbp), %eax
	cmpl	-16(%rbp), %eax
	jne	LBB1_2
## %bb.1:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -8(%rbp)
	jmp	LBB1_25
LBB1_2:
	movq	_seg@GOTPCREL(%rip), %rax
	movq	32(%rax), %rax
	movslq	-12(%rbp), %rcx
	movss	(%rax,%rcx,4), %xmm0    ## xmm0 = mem[0],zero,zero,zero
	cvtss2sd	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	jne	LBB1_3
	jp	LBB1_3
	jmp	LBB1_4
LBB1_3:
	movq	_seg@GOTPCREL(%rip), %rax
	movq	32(%rax), %rax
	movslq	-16(%rbp), %rcx
	movss	(%rax,%rcx,4), %xmm0    ## xmm0 = mem[0],zero,zero,zero
	cvtss2sd	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	jne	LBB1_5
	jp	LBB1_5
LBB1_4:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -8(%rbp)
	jmp	LBB1_25
LBB1_5:
	movq	_MinOvlpDist@GOTPCREL(%rip), %rax
	movq	_seg@GOTPCREL(%rip), %rcx
	movsd	LCPI1_0(%rip), %xmm0    ## xmm0 = mem[0],zero
	movq	8(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movaps	%xmm0, %xmm3
	mulsd	%xmm2, %xmm3
	addsd	%xmm3, %xmm1
	movsd	%xmm1, -40(%rbp)
	movq	8(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movaps	%xmm0, %xmm3
	mulsd	%xmm2, %xmm3
	subsd	%xmm3, %xmm1
	movsd	%xmm1, -48(%rbp)
	movq	8(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movaps	%xmm0, %xmm3
	mulsd	%xmm2, %xmm3
	addsd	%xmm3, %xmm1
	movsd	%xmm1, -56(%rbp)
	movq	8(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	mulsd	%xmm2, %xmm0
	subsd	%xmm0, %xmm1
	movsd	%xmm1, -64(%rbp)
	movq	8(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	8(%rdx), %xmm0          ## xmm0 = mem[0],zero,zero,zero
	movq	8(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	subss	8(%rdx), %xmm0
	movq	8(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	8(%rdx), %xmm1          ## xmm1 = mem[0],zero,zero,zero
	movq	8(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	subss	8(%rdx), %xmm1
	mulss	%xmm1, %xmm0
	movq	8(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	4(%rdx), %xmm1          ## xmm1 = mem[0],zero,zero,zero
	movq	8(%rcx), %rdx
	movslq	-16(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	subss	4(%rdx), %xmm1
	movq	8(%rcx), %rdx
	movslq	-12(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	4(%rdx), %xmm2          ## xmm2 = mem[0],zero,zero,zero
	movq	8(%rcx), %rcx
	movslq	-16(%rbp), %rdx
	imulq	$12, %rdx, %rdx
	addq	%rdx, %rcx
	subss	4(%rcx), %xmm2
	mulss	%xmm2, %xmm1
	addss	%xmm1, %xmm0
	cvtss2sd	%xmm0, %xmm0
	movsd	%xmm0, -32(%rbp)
	movsd	-32(%rbp), %xmm0        ## xmm0 = mem[0],zero
	cvtsd2ss	%xmm0, %xmm0
	sqrtss	%xmm0, %xmm0
	ucomiss	(%rax), %xmm0
	ja	LBB1_20
## %bb.6:
	movss	LCPI1_1(%rip), %xmm0    ## xmm0 = mem[0],zero,zero,zero
	movq	_MT@GOTPCREL(%rip), %rax
	movsd	-32(%rbp), %xmm1        ## xmm1 = mem[0],zero
	cvtsd2ss	%xmm1, %xmm1
	sqrtss	%xmm1, %xmm1
	mulss	192(%rax), %xmm0
	ucomiss	%xmm1, %xmm0
	jbe	LBB1_8
## %bb.7:
	movl	-12(%rbp), %edi
	callq	_segindex
	movl	-16(%rbp), %edi
	movl	%eax, -92(%rbp)         ## 4-byte Spill
	callq	_segindex
	movl	-92(%rbp), %ecx         ## 4-byte Reload
	cmpl	%eax, %ecx
	jne	LBB1_20
LBB1_8:
	movq	_bound@GOTPCREL(%rip), %rax
	cmpl	$1, (%rax)
	jne	LBB1_11
## %bb.9:
	movq	_lbound0@GOTPCREL(%rip), %rax
	movsd	-40(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm0, %xmm1
	ja	LBB1_20
## %bb.10:
	movq	_lbound0@GOTPCREL(%rip), %rax
	movsd	-56(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm0, %xmm1
	ja	LBB1_20
LBB1_11:
	movq	_bound@GOTPCREL(%rip), %rax
	cmpl	$1, 4(%rax)
	jne	LBB1_14
## %bb.12:
	movq	_rbound0@GOTPCREL(%rip), %rax
	movsd	-48(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	ja	LBB1_20
## %bb.13:
	movq	_rbound0@GOTPCREL(%rip), %rax
	movsd	-64(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	ja	LBB1_20
LBB1_14:
	movq	_bound@GOTPCREL(%rip), %rax
	cmpl	$6, 4(%rax)
	jne	LBB1_17
## %bb.15:
	movq	_rbound0@GOTPCREL(%rip), %rax
	movsd	-48(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	ja	LBB1_20
## %bb.16:
	movq	_rbound0@GOTPCREL(%rip), %rax
	movsd	-64(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	ja	LBB1_20
LBB1_17:
	movq	_bound@GOTPCREL(%rip), %rax
	cmpl	$6, (%rax)
	jne	LBB1_21
## %bb.18:
	movq	_lbound0@GOTPCREL(%rip), %rax
	movsd	-40(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm0, %xmm1
	ja	LBB1_20
## %bb.19:
	movq	_lbound0@GOTPCREL(%rip), %rax
	movsd	-56(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movss	(%rax), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	ucomisd	%xmm0, %xmm1
	jbe	LBB1_21
LBB1_20:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -24(%rbp)
	jmp	LBB1_22
LBB1_21:
	movsd	-40(%rbp), %xmm0        ## xmm0 = mem[0],zero
	cvtsd2ss	%xmm0, %xmm0
	movsd	-56(%rbp), %xmm1        ## xmm1 = mem[0],zero
	cvtsd2ss	%xmm1, %xmm1
	callq	_maxf
	cvtss2sd	%xmm0, %xmm0
	movsd	%xmm0, -72(%rbp)
	movsd	-48(%rbp), %xmm0        ## xmm0 = mem[0],zero
	cvtsd2ss	%xmm0, %xmm0
	movsd	-64(%rbp), %xmm1        ## xmm1 = mem[0],zero
	cvtsd2ss	%xmm1, %xmm1
	callq	_minf
	movq	_seg@GOTPCREL(%rip), %rax
	cvtss2sd	%xmm0, %xmm0
	movsd	%xmm0, -80(%rbp)
	movq	32(%rax), %rcx
	movslq	-12(%rbp), %rdx
	movss	(%rcx,%rdx,4), %xmm0    ## xmm0 = mem[0],zero,zero,zero
	movq	32(%rax), %rax
	movslq	-16(%rbp), %rcx
	addss	(%rax,%rcx,4), %xmm0
	cvtss2sd	%xmm0, %xmm0
	movsd	-72(%rbp), %xmm1        ## xmm1 = mem[0],zero
	subsd	-80(%rbp), %xmm1
	subsd	%xmm1, %xmm0
	movsd	%xmm0, -24(%rbp)
LBB1_22:
	movsd	LCPI1_2(%rip), %xmm0    ## xmm0 = mem[0],zero
	movsd	-24(%rbp), %xmm1        ## xmm1 = mem[0],zero
	ucomisd	%xmm0, %xmm1
	jbe	LBB1_24
## %bb.23:
	movsd	-24(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movsd	%xmm0, -8(%rbp)
	jmp	LBB1_25
LBB1_24:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -8(%rbp)
LBB1_25:
	movsd	-8(%rbp), %xmm0         ## xmm0 = mem[0],zero
	addq	$96, %rsp
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_xoverlap               ## -- Begin function xoverlap
	.p2align	4, 0x90
_xoverlap:                              ## @xoverlap
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$80, %rsp
	movl	%edi, -12(%rbp)
	movl	%esi, -16(%rbp)
	movl	-12(%rbp), %eax
	cmpl	-16(%rbp), %eax
	jne	LBB2_2
## %bb.1:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -8(%rbp)
	jmp	LBB2_3
LBB2_2:
	movq	_seg@GOTPCREL(%rip), %rax
	movl	-12(%rbp), %ecx
	shll	$1, %ecx
	subl	$1, %ecx
	movl	%ecx, -60(%rbp)
	movl	-12(%rbp), %ecx
	shll	$1, %ecx
	movl	%ecx, -64(%rbp)
	movl	-16(%rbp), %ecx
	shll	$1, %ecx
	subl	$1, %ecx
	movl	%ecx, -68(%rbp)
	movl	-16(%rbp), %ecx
	shll	$1, %ecx
	movl	%ecx, -72(%rbp)
	movq	8(%rax), %rdx
	movslq	-60(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm0           ## xmm0 = mem[0],zero,zero,zero
	cvtss2sd	%xmm0, %xmm0
	movq	8(%rax), %rdx
	movslq	-68(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rax), %rdx
	movslq	-60(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movq	32(%rax), %rax
	movslq	-68(%rbp), %rdx
	movss	(%rax,%rdx,4), %xmm3    ## xmm3 = mem[0],zero,zero,zero
	cvtss2sd	%xmm3, %xmm3
	callq	_sxoverlap
	movq	_seg@GOTPCREL(%rip), %rax
	movsd	%xmm0, -32(%rbp)
	movq	8(%rax), %rdx
	movslq	-60(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm0           ## xmm0 = mem[0],zero,zero,zero
	cvtss2sd	%xmm0, %xmm0
	movq	8(%rax), %rdx
	movslq	-72(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rax), %rdx
	movslq	-60(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movq	32(%rax), %rax
	movslq	-72(%rbp), %rdx
	movss	(%rax,%rdx,4), %xmm3    ## xmm3 = mem[0],zero,zero,zero
	cvtss2sd	%xmm3, %xmm3
	callq	_sxoverlap
	movq	_seg@GOTPCREL(%rip), %rax
	movsd	%xmm0, -40(%rbp)
	movq	8(%rax), %rdx
	movslq	-64(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm0           ## xmm0 = mem[0],zero,zero,zero
	cvtss2sd	%xmm0, %xmm0
	movq	8(%rax), %rdx
	movslq	-68(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rax), %rdx
	movslq	-64(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movq	32(%rax), %rax
	movslq	-68(%rbp), %rdx
	movss	(%rax,%rdx,4), %xmm3    ## xmm3 = mem[0],zero,zero,zero
	cvtss2sd	%xmm3, %xmm3
	callq	_sxoverlap
	movq	_seg@GOTPCREL(%rip), %rax
	movsd	%xmm0, -48(%rbp)
	movq	8(%rax), %rdx
	movslq	-64(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm0           ## xmm0 = mem[0],zero,zero,zero
	cvtss2sd	%xmm0, %xmm0
	movq	8(%rax), %rdx
	movslq	-72(%rbp), %rsi
	imulq	$12, %rsi, %rsi
	addq	%rsi, %rdx
	movss	(%rdx), %xmm1           ## xmm1 = mem[0],zero,zero,zero
	cvtss2sd	%xmm1, %xmm1
	movq	32(%rax), %rdx
	movslq	-64(%rbp), %rsi
	movss	(%rdx,%rsi,4), %xmm2    ## xmm2 = mem[0],zero,zero,zero
	cvtss2sd	%xmm2, %xmm2
	movq	32(%rax), %rax
	movslq	-72(%rbp), %rdx
	movss	(%rax,%rdx,4), %xmm3    ## xmm3 = mem[0],zero,zero,zero
	cvtss2sd	%xmm3, %xmm3
	callq	_sxoverlap
	movsd	%xmm0, -56(%rbp)
	movsd	-32(%rbp), %xmm0        ## xmm0 = mem[0],zero
	addsd	-40(%rbp), %xmm0
	addsd	-48(%rbp), %xmm0
	addsd	-56(%rbp), %xmm0
	movsd	%xmm0, -24(%rbp)
	movsd	-24(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movsd	%xmm0, -8(%rbp)
LBB2_3:
	movsd	-8(%rbp), %xmm0         ## xmm0 = mem[0],zero
	addq	$80, %rsp
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.section	__TEXT,__literal8,8byte_literals
	.p2align	3               ## -- Begin function sxoverlap
LCPI3_0:
	.quad	4562254508917369340     ## double 0.001
LCPI3_1:
	.quad	4602678819172646912     ## double 0.5
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_sxoverlap
	.p2align	4, 0x90
_sxoverlap:                             ## @sxoverlap
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$96, %rsp
	movsd	LCPI3_1(%rip), %xmm4    ## xmm4 = mem[0],zero
	movsd	%xmm0, -16(%rbp)
	movsd	%xmm1, -24(%rbp)
	movsd	%xmm2, -32(%rbp)
	movsd	%xmm3, -40(%rbp)
	movsd	-16(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movaps	%xmm4, %xmm1
	mulsd	-32(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -56(%rbp)
	movsd	-16(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movaps	%xmm4, %xmm1
	mulsd	-32(%rbp), %xmm1
	subsd	%xmm1, %xmm0
	movsd	%xmm0, -64(%rbp)
	movsd	-24(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movaps	%xmm4, %xmm1
	mulsd	-40(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -72(%rbp)
	movsd	-24(%rbp), %xmm0        ## xmm0 = mem[0],zero
	mulsd	-40(%rbp), %xmm4
	subsd	%xmm4, %xmm0
	movsd	%xmm0, -80(%rbp)
	movsd	-56(%rbp), %xmm0        ## xmm0 = mem[0],zero
	cvtsd2ss	%xmm0, %xmm0
	movsd	-72(%rbp), %xmm1        ## xmm1 = mem[0],zero
	cvtsd2ss	%xmm1, %xmm1
	callq	_maxf
	cvtss2sd	%xmm0, %xmm0
	movsd	%xmm0, -88(%rbp)
	movsd	-64(%rbp), %xmm0        ## xmm0 = mem[0],zero
	cvtsd2ss	%xmm0, %xmm0
	movsd	-80(%rbp), %xmm1        ## xmm1 = mem[0],zero
	cvtsd2ss	%xmm1, %xmm1
	callq	_minf
	movsd	LCPI3_0(%rip), %xmm1    ## xmm1 = mem[0],zero
	cvtss2sd	%xmm0, %xmm0
	movsd	%xmm0, -96(%rbp)
	movsd	-32(%rbp), %xmm0        ## xmm0 = mem[0],zero
	addsd	-40(%rbp), %xmm0
	movsd	-88(%rbp), %xmm2        ## xmm2 = mem[0],zero
	subsd	-96(%rbp), %xmm2
	subsd	%xmm2, %xmm0
	movsd	%xmm0, -48(%rbp)
	movsd	-48(%rbp), %xmm0        ## xmm0 = mem[0],zero
	ucomisd	%xmm1, %xmm0
	jbe	LBB3_2
## %bb.1:
	movsd	-48(%rbp), %xmm0        ## xmm0 = mem[0],zero
	movsd	%xmm0, -8(%rbp)
	jmp	LBB3_3
LBB3_2:
	xorps	%xmm0, %xmm0
	movsd	%xmm0, -8(%rbp)
LBB3_3:
	movsd	-8(%rbp), %xmm0         ## xmm0 = mem[0],zero
	addq	$96, %rsp
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.comm	_ovlp,48,3              ## @ovlp
	.comm	_seg,48,3               ## @seg
	.comm	_MinOvlpDist,4,2        ## @MinOvlpDist
	.comm	_MT,200,3               ## @MT
	.comm	_bound,8,2              ## @bound
	.comm	_lbound0,4,2            ## @lbound0
	.comm	_rbound0,4,2            ## @rbound0
	.comm	_perc_prob,4,2          ## @perc_prob
	.comm	_degprob,4,2            ## @degprob
	.comm	_mxmts,4,2              ## @mxmts
	.comm	_mxmotors,4,2           ## @mxmotors
	.comm	_twopop,4,2             ## @twopop
	.comm	_popdist,4,2            ## @popdist
	.comm	_t_length,4,2           ## @t_length
	.comm	_manual,4,2             ## @manual
	.comm	_FixedDensity,4,2       ## @FixedDensity
	.comm	_lengthmix,4,2          ## @lengthmix
	.comm	_TotFluxConst,4,2       ## @TotFluxConst
	.comm	_MOTOR,4,2              ## @MOTOR
	.comm	_D3_output,4,2          ## @D3_output
	.comm	_ACTIN,4,2              ## @ACTIN
	.comm	_AvNeighboredMTs,4,2    ## @AvNeighboredMTs
	.comm	_addmodel,4,2           ## @addmodel
	.comm	_polaritymodel,4,2      ## @polaritymodel
	.comm	_dirname,8,3            ## @dirname
	.comm	_MTfirst,4,2            ## @MTfirst
	.comm	_newMTs0,4,2            ## @newMTs0
	.comm	_newfreq,4,2            ## @newfreq
	.comm	_sum,48,3               ## @sum
	.comm	_track,144,3            ## @track
	.comm	_corr_ti,8,3            ## @corr_ti
	.comm	_corr_tf,8,3            ## @corr_tf
	.comm	_Lnet,8,2               ## @Lnet
	.comm	_reflectb0,4,2          ## @reflectb0
	.comm	_wfreq,4,2              ## @wfreq
	.comm	_PolarityRatio0,4,2     ## @PolarityRatio0
	.comm	_PolarityRatioNew,4,2   ## @PolarityRatioNew
	.comm	_ProbBipolar,4,2        ## @ProbBipolar
	.comm	_ProbBundling,4,2       ## @ProbBundling
	.comm	_ProbLegDown,4,2        ## @ProbLegDown
	.comm	_ProbActive,4,2         ## @ProbActive
	.comm	_ProbKinesin,4,2        ## @ProbKinesin
	.comm	_pFlp,4,2               ## @pFlp
	.comm	_neqsteps,4,2           ## @neqsteps
	.comm	_DomainX,4,2            ## @DomainX
	.comm	_DomainZ,4,2            ## @DomainZ
	.comm	_fstall_k,4,2           ## @fstall_k
	.comm	_fstall_d,4,2           ## @fstall_d
	.comm	_fstall_bp,4,2          ## @fstall_bp
	.comm	_vel0_k,4,2             ## @vel0_k
	.comm	_vel0_d,4,2             ## @vel0_d
	.comm	_vel0_bp,4,2            ## @vel0_bp
	.comm	_xim0,4,2               ## @xim0
	.comm	_xisol,4,2              ## @xisol
	.comm	_sticky_drag,4,2        ## @sticky_drag
	.comm	_lamb_dyn,4,2           ## @lamb_dyn
	.comm	_dist_dyn,8,3           ## @dist_dyn
	.comm	_lamb_kin,4,2           ## @lamb_kin
	.comm	_dist_kin,8,3           ## @dist_kin
	.comm	_ProbBi_dist,8,3        ## @ProbBi_dist
	.comm	_ProbKi_dist,8,3        ## @ProbKi_dist
	.comm	_ProbDy_dist,8,3        ## @ProbDy_dist
	.comm	_ProbAct_dist,8,3       ## @ProbAct_dist
	.comm	_vel_flux,4,2           ## @vel_flux
	.comm	_iter,4,2               ## @iter
	.comm	_niter,4,2              ## @niter
	.comm	_iav,4,2                ## @iav
	.comm	_nav,4,2                ## @nav
	.comm	_lboundNew,4,2          ## @lboundNew
	.comm	_rboundNew,4,2          ## @rboundNew
	.comm	_dt,4,2                 ## @dt
	.comm	_exclude,4,2            ## @exclude
	.comm	_nbox,4,2               ## @nbox
	.comm	_nbox2,4,2              ## @nbox2
	.comm	_hist_vel_min,4,2       ## @hist_vel_min
	.comm	_hist_vel_max,4,2       ## @hist_vel_max
	.comm	_hist_absvel_min,4,2    ## @hist_absvel_min
	.comm	_hist_absvel_max,4,2    ## @hist_absvel_max
	.comm	_randav,4,2             ## @randav
	.comm	_randstd,4,2            ## @randstd
	.comm	_NmaxMTs,4,2            ## @NmaxMTs
	.comm	_velpol,4,2             ## @velpol
	.comm	_veldepol,4,2           ## @veldepol
	.comm	_ext_fr,4,2             ## @ext_fr
	.comm	_ext_fl,4,2             ## @ext_fl
	.comm	_touch_depth,4,2        ## @touch_depth
	.comm	_spring_r,4,2           ## @spring_r
	.comm	_spring_l,4,2           ## @spring_l
	.comm	_force_it0,4,2          ## @force_it0
	.comm	_actin,8,2              ## @actin
	.comm	_av,648,3               ## @av
	.comm	_bundle,20,2            ## @bundle
	.comm	_maxvel,4,2             ## @maxvel
	.comm	_idum,8,3               ## @idum
	.comm	_idum2,8,3              ## @idum2
	.comm	_kBT,4,2                ## @kBT
	.comm	_lin_density,4,2        ## @lin_density
	.comm	_forcepol,4,2           ## @forcepol
	.comm	_forcelength,4,2        ## @forcelength
	.comm	_MCfreq,4,2             ## @MCfreq
.subsections_via_symbols
