✈️ TweetUS Airline – Analisi Sentiment con Spark

📚 Descrizione del Progetto

Questo progetto analizza 14.640 tweet pubblicati negli Stati Uniti nel febbraio 2015, relativi a compagnie aeree americane. L'obiettivo è esplorare, visualizzare e modellare i dati tramite **Apache Spark**, utilizzando sia **Spark SQL** per l'analisi esplorativa sia **Spark ML** per la modellazione predittiva.


🛠️ Tecnologie Utilizzate

- **Apache Spark** (SQL, MLlib)
- **Python / PySpark**
- **Librerie esterne**: `matplotlib`, `wordcloud`, `Pandas`
- Dataset in formato **CSV**



🧩 Struttura del Dataset

- **Totale osservazioni**: 14.640 tweet
- **Variabili principali**:
  - `airline_sentiment`: sentiment (positive, neutral, negative)
  - `airline`: compagnia aerea
  - `name`: autore del tweet
  - `retweet_count`: numero di retweet
  - `text`: contenuto testuale
  - `negativereason`: motivazione del sentiment negativo
  - `user_timezone`: fuso orario utente



🔍 Analisi Esplorativa (Spark SQL)

1. Tweet per utente
- Raggruppamento per `name`, ordinamento decrescente e visualizzazione dei top 10 utenti.

2. Retweet per compagnia
- Due approcci:
  - Filtraggio su `retweet_count > 0` prima dell’aggregazione
  - Aggregazione con `sum` su `retweet_count`, filtrando successivamente

3. Tag Cloud dei commenti
- WordCloud generata con libreria `wordcloud` da tutti i `text` presenti.

4. Compagnia con più tweet negativi
- Filtraggio su `airline_sentiment == 'negative'` e conteggio per `airline`

5. Distribuzione sentiment per compagnia
- Raggruppamento combinato `airline` + `airline_sentiment` con `count`

6. Distribuzione motivazioni negative
- Aggregazione su `negativereason` filtrando tweet con sentiment negativo

7. Autori con >10 retweet
- Somma di `retweet_count` per autore, filtrati su soglia >10 e ordinati alfabeticamente

8–9. WordCloud per sentiment negativo/positivo dato il nome compagnia



🤖 Analisi Predittiva (Spark ML)

📌 Obiettivo 1: **Classificazione del sentiment**

🔧 Tecniche usate:
- **Tokenizer** e **StopWordsRemover** su `text` e `negativereason`
- **CountVectorizer** + **IDF** per estrarre feature testuali
- **StringIndexer** per target `airline_sentiment` e `airline`
- **OneHotEncoder** per feature categoriche
- **VectorAssembler** per creare un'unica feature vector
- **Modello**: `LogisticRegression`
- **Valutazione**: Matrice di confusione, Precision, Recall, Accuracy

🎯 Classi previste:
- `positive`, `neutral`, `negative` (indicizzate con StringIndexer)



📌 Obiettivo 2: **Regressione sul numero di retweet**

- Stessa pipeline di preprocessing
- **Modello**: `LinearRegression`
- **Target**: `retweet_count`
- **Valutazione**: RMSE (Root Mean Squared Error), MAE (Mean Absolute Error)


🎁 Bonus: Classificazione con n-gram e Naive Bayes

- Introduzione del **n-gram** (bigrammi) nella pipeline di trasformazione
- Nuovo classificatore: `NaiveBayes` (multinomial)
- Pipeline aggiornata con output WordCloud, matrice di confusione e metriche di performance



📊 Risultati finali

- Gli utenti più attivi si trovano nella zona **Eastern Time (US & Canada)**
- La **compagnia United** riceve il maggior numero di tweet negativi
- La **motivazione negativa più comune** è "Customer Service Issue"
- Il **sentiment è prevedibile con buoni risultati** tramite Logistic Regression
- Il **numero di retweet** può essere stimato con discreta accuratezza

