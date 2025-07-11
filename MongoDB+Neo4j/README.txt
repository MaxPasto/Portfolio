# TweetUS Airline - Analisi Sentiment su Compagnie Aeree USA

## 📌 Descrizione del Progetto

Il progetto analizza 14.848 tweet raccolti nel febbraio 2015 relativi a compagnie aeree statunitensi. L'analisi è stata condotta sia su **MongoDB (NoSQL orientato ai documenti)** che su **Neo4j (database a grafo)**, mettendo in luce diverse prospettive attraverso interrogazioni aggregate e visualizzazioni.

---

## 📦 Dataset

- **Periodo**: Febbraio 2015  
- **Localizzazione geografica**: USA  
- **Totale tweet**: 14.848  
- **Formato origine**: CSV

---

## 🧩 Struttura MongoDB

- **Database**: `Tweet`
- **Collezione**: `tweet`

### Variabili principali utilizzate
- `airline_sentiment`: classificazione del sentiment (positive, neutral, negative)
- `airline`: compagnia aerea
- `name`: nome dell’utente
- `retweet_count`: numero di retweet
- `text`: contenuto del tweet
- `negativereason`: motivazione del sentiment negativo
- `user_timezone`: fuso orario dell’utente

### Analisi MongoDB implementate

1. 🔢 Numero di tweet per utente
2. 🔁 Numero di retweet per compagnia aerea
3. ☁️ WordCloud dei termini più usati nei tweet
4. 😡 Compagnia con più tweet negativi
5. 📊 Distribuzione del sentiment per compagnia
6. ❗ Distribuzione delle motivazioni negative
7. ✍️ Utenti con più di 10 retweet ricevuti
8. ☁️ WordCloud per sentiment negativo dato il nome compagnia
9. ☁️ WordCloud per sentiment positivo dato il nome compagnia

---

## 🌐 Struttura Neo4j

- **Entità (nodi)**:
  - `Tweet`
  - `Users`
  - `Airlines`
- **Relazioni**:
  - `HASCRITTO`: un utente ha scritto un tweet
  - `RETWEETED`: un tweet è stato ritwittato
  - `RIFERITO`: il tweet si riferisce a una compagnia aerea

### Query principali in Cypher

1. 👤 Numero di tweet per utente
2. 🔁 Totale retweet per compagnia
3. ❌ Compagnia con più tweet negativi
4. 📊 Distribuzione sentiment per compagnia
5. ❗ Distribuzione motivazioni negative
6. ✍️ Autori con più di 10 retweet ricevuti

---

## 🆕 Estensioni aggiuntive (Neo4j)

1. 🌎 Distribuzione dei tweet negativi per compagnia e fuso orario
2. 🕒 Analisi delle motivazioni negative per fuso orario "Eastern Time (US & Canada)"
3. 🛫 Analisi delle motivazioni negative per la compagnia **United**

---

## 🛠️ Librerie e strumenti utilizzati

- **MongoDB** e **MongoDB Compass**
- **Neo4j**
- **Python**:
  - `pymongo`
  - `matplotlib`
  - `wordcloud`
  - `Pandas`, `NumPy`
- **Cypher** per le query grafiche in Neo4j

---

## 📌 Note finali

Il progetto mette in evidenza:
- Il valore delle analisi esplorative sui dati social.
- L’utilizzo efficace di due approcci NoSQL complementari (documenti e grafi).
- L’impatto che un buon modello dati e query ben progettate hanno sulla capacità analitica.

---
