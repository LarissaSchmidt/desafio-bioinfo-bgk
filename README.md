# Desafio Bioinformática - BGK
![WhatsApp-Image-2020-12-03-at-11 25 31](https://github.com/user-attachments/assets/349ed1d1-1749-4723-93b0-dfd6dd9704d8)
Desafio Técnico de BioInformática proposto pela empresa Biogenetika como parte do processo de entrevista para emprego. 
# Badges 
![Badge em Desenvolvimento]( https://img.shields.io/badge/Em%20Análise%20-8A2BE2)

# Descrição do Projeto

Está é a resolução de um desafio Técnico de BioInformática proposto pela empresa Biogenetika como parte do processo de entrevista para emprego. 
O desafio foi feito em linguagem Python e baseado em um problema retirado do site Rosalind, onde deveria ser realizada uma recuperação de dados de dentro do site NCBI.

O problema pode ser acessado no link: <https://rosalind.info/problems/frmt/>

Nesse mesmo programa, houve a tentativa de automatizar a entrada de dados como um parâmetro que seria colocado pelo próprio usuário no terminal utilizado.  Também foi colocado linhas para tratamento básico diante de erros que poderiam surgir, como a inserção de uma ID inválida. 

# Construção código

Pensando no objetivo do código, que é: acessar o banco de dados, pegar uma sequência de nucleotídeos referentes a uma ID específica, ler o arquivo FASTA e retornar a sequência referente. Prosseguimos com os seguintes passos:   

Definimos uma função com 'def', para facilitar a visualização e para utilizar posteriormente no nosso código e outra 'genbank_id' para identificar as futuras IDs
```ruby
def pega_seq_genbank(genbank_id):
```

Para acessar o banco de dados de nucleotideos do NCBI no site: https://www.ncbi.nlm.nih.gov/nucleotide/

Precisamos instalar a library Biopython dentro do console e importar duas funcionalidades:
```ruby
from Bio import Entrez, SeqIO
```

A função Entrez acessa o site pelo www, nos permitindo trabalhar com o conjunto de variáveis necessárias e acessar os url das IDs. O acesso se da pelo login de e-mail setado (Já previamente cadastrado no NCBI)
```ruby
Entrez.email = 'schmidttein@gmail.com'
```

A função efetch é a que retorna do banco de dados as informações e nos permite lidar com as ferramentas como um objeto 
```ruby
Entrez.efetch
```

Utilizamos a função SeqIO.read para ler o arquivo em formato FASTA das sequências nucleotídicas contidas na variável 'conector' e suas respectivas IDs
```ruby
SeqIO.read(conector, 'fasta')
```

Juntando essas funções e definindo as variáveis de 'registros', 'conector' e 'entrada' temos:
```ruby
def pega_seq_genbank(genbank_id):
        from Bio import Entrez, SeqIO
        Entrez.email = 'schmidttein@gmail.com'
        registros = []
        conector = Entrez.efetch(db='nucleotide', id=genbank_id, rettype='fasta')
        entrada = SeqIO.read(conector, 'fasta')
        registros.append(entrada)
        conector.close()
```
A linha com 'conector.close()' finaliza a leitura do arquivo.

Com isso acessamos o banco de dados e suas informações, para padronizar o tamanho e pegar a menor sequência referente a determinada ID
Definimos então uma variável 'menor_registro' que servirá para definir a menor sequência nucleotídica referente aquela ID. É definida pelo minímo 'min', lenght 'len' dentro de 'registros'.

```ruby
menor_registro = min(registros, key=lambda x: len(x.seq))
 tamanho = len(menor_registro.seq)
```

Para caráter de exemplificação, optei por organizar as informações da sequência obtida em uma string única contida na variável 'seq_fasta'

```ruby
        seq_fasta = (f'O código >{menor_registro.id} é referente a uma sequência com {tamanho} pares de base '
                     f'tal sequência é: '
                     f'\n{menor_registro.seq}')
        return seq_fasta
```

Até este momento definimos 'def' a função 'pega_seq_genbank' que utilizaremos para realizar a tarefa
Agora perguntaremos para o usuário de qual ID ele gostaria de obter a sequência FASTA. Para tal, utilizamos um 'input' que permite o usuário digitar a ID desejada
```ruby
genbank_id = [str(input('Digite o id do genbank: ').upper().strip())]
```
Em seguida, definimos uma variável 'final' para guardar a string a ser printada. Chamamos a função 'pega_seq_genbank', antes definida e a rodamos com a variável 'genbank_id' colocada pelo usuário e a printamos na tela a variável 'final' com a função 'print'
```ruby
genbank_id = [str(input('Digite o id do genbank: ').upper().strip())]
    final = pega_seq_genbank(genbank_id)
    print(final)
```
Utilizamos as funções 'while True' e 'while _ not in _' para que o programa possa rodar quantas vezes o usuário pedir. A variável 'sn' é utilizada somente para checar a escolha do usuário de continuar ou não; e caso não queira dar sequência ao comando 'break' assim encerrando o programa com o 'print('Até mais!')'

```ruby
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
```

Para evitar os erros que poderiam advir de uma ID inexistente ou errada, realizamos um tratamento de erros com o comando 'try' e 'except'; que tenta rodar o código e caso de erro retorna a mensagem de erro sem finalizar o programa, continuando o loop normalmente para receber o próximo input de ID
```ruby
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
```

Com isso, temos o nosso código completo

# Código

```ruby
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

```

# Índice 

* [Badges](#badges)
* [Índice](#índice)
* [Descrição do Projeto](#descrição-do-projeto)
* [Construção Código](#construção-código)
* [Código](#código)
* [Explicação](#explicação)
* [Tecnologias utilizadas](#tecnologias-utilizadas)
* [Pessoas Contribuidoras](#pessoas-contribuidoras)
* [Pessoas Desenvolvedoras do Projeto](#pessoas-desenvolvedoras)
* [Licença](#licença)
* [Conclusão](#conclusão)
