def pega_seq_genbank(genbank_id):
    try:
        from Bio import Entrez, SeqIO
        Entrez.email = 'schmidttein@gmail.com'
        registros = []
        conector = Entrez.efetch(db='nucleotide', id=genbank_id, rettype='fasta')
        entrada = SeqIO.read(conector, 'fasta')
        registros.append(entrada)
        conector.close()
        menor_registro = min(registros, key=lambda x: len(x.seq))
        tamanho = len(menor_registro.seq)
        seq_fasta = (f'O código >{menor_registro.id} é referente a uma sequência com {tamanho} pares de base '
                     f'tal sequência é: '
                     f'\n{menor_registro.seq}')
        return seq_fasta
    except Exception as c:
        print(f"Erro ao buscar sequência: {c}")
        return None


while True:
    genbank_id = [str(input('Digite o id do genbank: ').upper().strip())]
    final = pega_seq_genbank(genbank_id)
    print(final)
    sn = ' '
    while sn not in 'SN':
        sn = input('Deseja ler outra sequencia: [S/N] ').strip().upper()[0]
    if sn == 'N':
        break
print('Até mais!')
