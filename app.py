import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="Protein Structure Predictor", page_icon="🧬", layout="centered")

# Chou-Fasman propensity values: (Alpha Helix, Beta Sheet)
propensity = {
    "A": (1.42, 0.83), "R": (0.98, 0.93), "N": (0.67, 0.89),
    "D": (1.01, 0.54), "C": (0.70, 1.19), "Q": (1.11, 1.10),
    "E": (1.51, 0.37), "G": (0.57, 0.75), "H": (1.00, 0.87),
    "I": (1.08, 1.60), "L": (1.21, 1.30), "K": (1.14, 0.74),
    "M": (1.45, 1.05), "F": (1.13, 1.38), "P": (0.57, 0.55),
    "S": (0.77, 0.75), "T": (0.83, 1.19), "W": (1.08, 1.37),
    "Y": (0.69, 1.47), "V": (1.06, 1.70)
}

def validate_sequence(sequence):
    valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return bool(sequence) and all(amino_acid in valid_amino_acids for amino_acid in sequence)

def predict_structure(sequence):
    n = len(sequence)
    helix_region = [False] * n
    sheet_region = [False] * n

    # STEP 1: Find alpha helix nucleation regions
    for i in range(n - 5):
        window = sequence[i:i+6]
        count = sum(1 for amino_acid in window if propensity[amino_acid][0] > 1.0)

        if count >= 4:
            start, end = i, i + 6

            while start > 0:
                values = [propensity[a][0] for a in sequence[start-1:start+3]]
                if sum(values) / 4 >= 1.0:
                    start -= 1
                else:
                    break

            while end < n:
                values = [propensity[a][0] for a in sequence[end-3:end+1]]
                if sum(values) / 4 >= 1.0:
                    end += 1
                else:
                    break

            region = sequence[start:end]
            alpha_avg = sum(propensity[a][0] for a in region) / len(region)
            beta_avg = sum(propensity[a][1] for a in region) / len(region)

            if len(region) >= 6 and alpha_avg > 1.03 and alpha_avg > beta_avg:
                for j in range(start, end):
                    helix_region[j] = True

    # STEP 2: Find beta sheet nucleation regions
    for i in range(n - 4):
        window = sequence[i:i+5]
        count = sum(1 for amino_acid in window if propensity[amino_acid][1] > 1.0)

        if count >= 3:
            start, end = i, i + 5

            while start > 0:
                values = [propensity[a][1] for a in sequence[start-1:start+3]]
                if sum(values) / 4 >= 1.0:
                    start -= 1
                else:
                    break

            while end < n:
                values = [propensity[a][1] for a in sequence[end-3:end+1]]
                if sum(values) / 4 >= 1.0:
                    end += 1
                else:
                    break

            region = sequence[start:end]
            alpha_avg = sum(propensity[a][0] for a in region) / len(region)
            beta_avg = sum(propensity[a][1] for a in region) / len(region)

            if len(region) >= 5 and beta_avg > 1.05 and beta_avg > alpha_avg:
                for j in range(start, end):
                    sheet_region[j] = True

    # STEP 3: Generate final prediction
    prediction = ""
    for i in range(n):
        if helix_region[i] and sheet_region[i]:
            prediction += "H" if propensity[sequence[i]][0] > propensity[sequence[i]][1] else "E"
        elif helix_region[i]:
            prediction += "H"
        elif sheet_region[i]:
            prediction += "E"
        else:
            prediction += "C"

    return prediction

st.title("🧬 Protein Secondary Structure Predictor")
st.write("Predict Alpha Helix, Beta Sheet, and Coil regions using the Chou–Fasman method.")
st.info("Enter a protein sequence using the 20 standard one-letter amino acid codes.")

sequence_input = st.text_area(
    "Protein sequence",
    placeholder="Example: MKTFFVAGL...",
    height=120
)

if st.button("Predict structure", type="primary"):
    sequence = "".join(sequence_input.split()).upper()

    if not validate_sequence(sequence):
        st.error("Invalid sequence. Use only the 20 standard amino acid letters: ACDEFGHIKLMNPQRSTVWY.")
    else:
        result = predict_structure(sequence)
        total = len(result)
        h_count = result.count("H")
        e_count = result.count("E")
        c_count = result.count("C")

        st.success("Prediction completed!")
        st.subheader("Prediction result")
        st.write("**Input sequence**")
        st.code(sequence)
        st.write("**Predicted structure**")
        st.code(result)

        st.subheader("Prediction summary")
        col1, col2, col3 = st.columns(3)
        col1.metric("Alpha Helix (H)", h_count, f"{h_count / total * 100:.2f}%")
        col2.metric("Beta Sheet (E)", e_count, f"{e_count / total * 100:.2f}%")
        col3.metric("Coil (C)", c_count, f"{c_count / total * 100:.2f}%")

        fig, ax = plt.subplots()
        structures = ["Alpha Helix", "Beta Sheet", "Coil"]
        counts = [h_count, e_count, c_count]
        ax.bar(structures, counts)
        ax.set_title("Protein Secondary Structure Prediction")
        ax.set_xlabel("Secondary Structure")
        ax.set_ylabel("Number of Amino Acids")
        st.pyplot(fig)
        plt.close(fig)
