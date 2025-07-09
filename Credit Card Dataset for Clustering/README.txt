Analisi Cluster su Dataset di Carte di Credito


Descrizione del progetto

Questo progetto ha come obiettivo l'analisi di un dataset relativo all'utilizzo delle carte di credito, con particolare attenzione al comportamento dei clienti. 

L'analisi è incentrata su tecniche di clustering per identificare gruppi di clienti con caratteristiche simili.

Tecniche utilizzate

- K-Means Clustering  
  Metodo di clustering basato sulla minimizzazione della distanza intra-cluster. È stato utilizzato per individuare gruppi omogenei di clienti sulla base delle variabili principali del dataset.

- K-Medoids Clustering  
  Variante più robusta del K-Means, che seleziona osservazioni reali come centroidi (medoids). È particolarmente utile in presenza di outlier.

- Clustering Gerarchico Agglomerativo  
  Tecnica basata sulla fusione successiva dei cluster più simili. È stato costruito un dendrogramma per esplorare visivamente la struttura gerarchica dei dati e scegliere il numero di cluster.

- PCA (Principal Component Analysis)
  Tecnica di riduzione della dimensionalità usata per trasformare le variabili originali in componenti principali non correlate.  
  L’analisi è stata condotta anche nello spazio ridotto delle prime due componenti principali, per facilitare la visualizzazione dei cluster e per migliorare la qualità del raggruppamento.

Obiettivo

Raggruppare i clienti in cluster significativi sulla base dei loro comportamenti di acquisto, saldo, frequenza e uso del credito. 
Questo tipo di segmentazione può essere utile, ad esempio, per strategie di marketing mirate o per la gestione del rischio di credito.

Visualizzazioni

Sono stati prodotti grafici per:

- Valutare la bontà del clustering tramite Silhouette Score e metodo del gomito.
- Visualizzare i cluster nello spazio delle componenti principali (PC1 e PC2).
- Analizzare la struttura gerarchica tramite dendrogrammi.