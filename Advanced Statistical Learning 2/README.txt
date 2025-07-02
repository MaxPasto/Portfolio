Questo repository contiene due progetti sviluppati nell'ambito del corso Advanced Statistical Learning 2, focalizzati sull'analisi e la
previsione mediante reti neurali. 

I due casi studio affrontano problematiche reali su dati economici e commerciali, sfruttando tecniche di apprendimento statistico avanzato.

🔨 Caso Studio 1: Indice di Produzione delle Costruzioni – Italia e Francia
📄 Visualizza il progetto online: 

https://maxpasto.github.io/Portfolio/Advanced%20Statistical%20Learning%202/Indice%20di%20produzione%20delle%20costruzioni%20per%20Italia%20e%20Francia.html

👉 Indice di produzione delle costruzioni – HTML

📊 Fonte dei dati: https://ec.europa.eu/eurostat/databrowser/view/ei_isbu_m/default/table?lang=en
👉 Eurostat – Indice produzione costruzioni

Descrizione:
Il progetto analizza le serie storiche mensili (gennaio 1995 – dicembre 2022) dell'indice di produzione del settore costruzioni per Italia e Francia. 
Le serie, destagionalizzate e normalizzate al livello 2015=100, sono utilizzate per:

1) Analisi descrittiva delle caratteristiche principali delle serie.

2) Implementazione di reti neurali per la previsione dell’andamento futuro.

Le osservazioni tra gennaio 1995 e dicembre 2019 sono utilizzate per il training, mentre quelle successive costituiscono il test set per valutare le 
capacità previsive del modello.

🚗 Caso Studio 2: Previsione del prezzo delle auto usate
📄 Visualizza il progetto online:

 https://maxpasto.github.io/Portfolio/Advanced%20Statistical%20Learning%202/Vehicle%20dataset.html
👉 Vehicle dataset – HTML

📊 Fonte dei dati: https://www.kaggle.com/datasets/nehalbirla/vehicle-dataset-from-cardekho?select=Car+details+v3.csv
👉 Kaggle – Vehicle Dataset

Descrizione:
Questo studio ha l'obiettivo di prevedere il prezzo di vendita di auto usate in base a variabili numeriche e categoriali, sfruttando:

1) Una rete neurale feedforward implementata con la libreria nnet.

2) Un confronto con un modello benchmark: la regressione lineare.

Il focus è valutare se la rete neurale riesca a cogliere meglio le relazioni non lineari tra la variabile dipendente (prezzo) e le covariate. 
Il dataset è stato pulito da valori mancanti, righe duplicate e variabili mal formattate per garantire coerenza e qualità nella modellazione.