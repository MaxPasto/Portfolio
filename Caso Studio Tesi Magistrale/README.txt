Volatilità Basata sul Range: Un'Analisi Comparativa
Descrizione del progetto
Caso studio volto a confrontare modelli previsivi di volatilità finanziaria 
su tre delle principali società incluse nell’indice S&P 500: Apple (AAPL), Amazon (AMZN) e NVIDIA (NVDA).
I dati sono stati raccolti da Yahoo Finance e coprono il periodo dal 1° gennaio 2021 al 31 dicembre 2024.

L'obiettivo dell'analisi è valutare se l'uso di misure di volatilità basate sul range, come la differenza tra il prezzo massimo e minimo giornaliero,
possa migliorare l'accuratezza previsionale rispetto ai modelli tradizionali basati sui soli rendimenti (open-to-close, close-to-close).

⚙️ Modelli implementati
Basati sui rendimenti (Close-to-Close)

GARCH(1,1)

EGARCH(1,1)

Basati sul range

RANGE-GARCH(1,1)

HAR (Heterogeneous Autoregressive Model)

NN-HAR: estensione non lineare del modello HAR con architettura neurale

🧪 Valutazione delle performance
Le performance previsive dei modelli vengono valutate su campioni out-of-sample di:

50 giorni (breve termine)

100 giorni (medio termine)

200 giorni (lungo termine)

La volatilità di riferimento ("vera") è stimata con lo stimatore di Parkinson, e l’errore di previsione viene misurato tramite RMSE (Root Mean Squared Error).

Per confronti statistici tra modelli, viene applicato il test di Diebold-Mariano, per verificare se le differenze di accuratezza previsiva siano significative.

🔍 Estensione: modello non lineare
Nella seconda parte dello studio viene proposta un'estensione non lineare del modello HAR, in cui si considera una possibile relazione 
non lineare tra la volatilità corrente e la volatilità laggata 
𝑉(𝑡−1), la media su 5 giorni (settimanale), V_t(5), la media su 22 giorni (mensile) 𝑉_t(22).
​
 L’obiettivo è esplorare se una struttura più flessibile possa migliorare ulteriormente le performance previsive.
