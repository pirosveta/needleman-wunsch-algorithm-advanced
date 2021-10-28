package ru.bmstu.bioinformatics;

public class PairAlignment {
    private final String OUTPUT_SEQUENCE_1 = "Seq1: ", OUTPUT_SEQUENCE_2 = "Seq2: ", OUTPUT_SCORE = "Score: ",
            NEXT_LINE = "\n";
    private final char GAP = '_';
    private final int ZERO = 0, GAP_INDEX = 0, ONE = 1, FIRST_INDEX = 1, NUMBER_OF_SYMBOLS_IN_LINE = 50,
            OPEN, INFINITY, EXTEND = -1;

    private Cell[][] matrixMatch, matrixInsertion, matrixDeletion;
    private String firstSequence,
            secondSequence;
    private ScoringFunction scoringFunction;

    private StringBuilder firstAlignSequence = new StringBuilder(),
            secondAlignSequence = new StringBuilder();
    private int score;

    private int getIndel() {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getIndel();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getIndel();
        } else return ((BLOSUM62) scoringFunction).getIndel();
    }

    private int getMatch(char firstChar, char secondChar) {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getMatch();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getMatch();
        } else return ((BLOSUM62) scoringFunction).getScore(firstChar, secondChar);
    }

    private int getMismatch(char firstChar, char secondChar) {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getMismatch();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getMismatch();
        } else return ((BLOSUM62) scoringFunction).getScore(firstChar, secondChar);
    }

    private int checkMatch(int lineIndex, int columnIndex) {
        char firstChar = firstSequence.charAt(columnIndex - ONE),
                secondChar = secondSequence.charAt(lineIndex - ONE);
        if (firstChar == secondChar) {
            return getMatch(firstChar, secondChar);
        } else return getMismatch(firstChar, secondChar);
    }

    private void fillGapCells() {
        matrixMatch[GAP_INDEX][GAP_INDEX] = new Cell(null, PredecessorType.NULL, GAP, GAP, ZERO);
        matrixInsertion[GAP_INDEX][GAP_INDEX] = new Cell(null, PredecessorType.NULL, GAP, GAP, INFINITY);
        matrixDeletion[GAP_INDEX][GAP_INDEX] = new Cell(null, PredecessorType.NULL, GAP, GAP, INFINITY);

        for (int columnIndex = FIRST_INDEX; columnIndex <= firstSequence.length(); columnIndex++) {
            matrixMatch[GAP_INDEX][columnIndex] = new Cell(
                    matrixMatch[GAP_INDEX][columnIndex - ONE],
                    PredecessorType.LEFT,
                    firstSequence.charAt(columnIndex - ONE),
                    GAP,
                    INFINITY);
            matrixInsertion[GAP_INDEX][columnIndex] = new Cell(
                    matrixInsertion[GAP_INDEX][columnIndex - ONE],
                    PredecessorType.LEFT,
                    firstSequence.charAt(columnIndex - ONE),
                    GAP,
                    OPEN + (columnIndex - ONE) * EXTEND);
            matrixDeletion[GAP_INDEX][columnIndex] = new Cell(
                    matrixDeletion[GAP_INDEX][columnIndex - ONE],
                    PredecessorType.LEFT,
                    firstSequence.charAt(columnIndex - ONE),
                    GAP,
                    INFINITY);
        }
        for (int lineIndex = FIRST_INDEX; lineIndex <= secondSequence.length(); lineIndex++) {
            matrixMatch[lineIndex][GAP_INDEX] = new Cell(
                    matrixMatch[lineIndex - ONE][GAP_INDEX],
                    PredecessorType.UP,
                    GAP,
                    secondSequence.charAt(lineIndex - ONE),
                    INFINITY);
            matrixInsertion[lineIndex][GAP_INDEX] = new Cell(
                    matrixInsertion[lineIndex - ONE][GAP_INDEX],
                    PredecessorType.UP,
                    GAP,
                    secondSequence.charAt(lineIndex - ONE),
                    INFINITY);
            matrixDeletion[lineIndex][GAP_INDEX] = new Cell(
                    matrixDeletion[lineIndex - ONE][GAP_INDEX],
                    PredecessorType.UP,
                    GAP,
                    secondSequence.charAt(lineIndex - ONE),
                    OPEN + (lineIndex - ONE) * EXTEND);
        }
    }

    private void fillMatchMatrix(int lineIndex, int columnIndex) {
        int valueMatch = matrixMatch[lineIndex - 1][columnIndex - 1].getValue(),
                valueInsertion = matrixInsertion[lineIndex - 1][columnIndex - 1].getValue(),
                valueDeletion = matrixDeletion[lineIndex - 1][columnIndex - 1].getValue();
        int valueMatchOrMismatch = checkMatch(lineIndex, columnIndex);
        MatrixType matrixType = MatrixType.MATCH;
        Cell predecessorCell;
        int predecessorValue;

        if (valueInsertion >= valueMatch && valueInsertion >= valueDeletion) {
            matrixType = MatrixType.INSERTION;
        } else if (valueDeletion >= valueMatch && valueDeletion >= valueInsertion) {
            matrixType = MatrixType.DELETION;
        }
        switch (matrixType) {
            case INSERTION:
                predecessorCell = matrixInsertion[lineIndex - ONE][columnIndex - ONE];
                predecessorValue = valueInsertion + valueMatchOrMismatch;
                break;
            case DELETION:
                predecessorCell = matrixDeletion[lineIndex - ONE][columnIndex - ONE];
                predecessorValue = valueDeletion + valueMatchOrMismatch;
                break;
            default:
                predecessorCell = matrixMatch[lineIndex - ONE][columnIndex - ONE];
                predecessorValue = valueMatch + valueMatchOrMismatch;
        }

        matrixMatch[lineIndex][columnIndex] = new Cell(
                predecessorCell,
                PredecessorType.DIAG,
                firstSequence.charAt(columnIndex - ONE),
                secondSequence.charAt(lineIndex - ONE),
                predecessorValue);
    }

    private void fillInsertionMatrix(int lineIndex, int columnIndex) {
        int valueMatch = matrixMatch[lineIndex][columnIndex - 1].getValue() + OPEN,
                valueInsertion = matrixInsertion[lineIndex][columnIndex - 1].getValue() + EXTEND,
                valueDeletion = matrixDeletion[lineIndex][columnIndex - 1].getValue() + OPEN;
        MatrixType matrixType = MatrixType.INSERTION;
        Cell predecessorCell;
        int predecessorValue;

        if (valueMatch >= valueInsertion && valueMatch >= valueDeletion) {
            matrixType = MatrixType.MATCH;
        } else if (valueDeletion >= valueMatch && valueDeletion >= valueInsertion) {
            matrixType = MatrixType.DELETION;
        }
        switch (matrixType) {
            case MATCH:
                predecessorCell = matrixMatch[lineIndex][columnIndex - ONE];
                predecessorValue = valueMatch;
                break;
            case DELETION:
                predecessorCell = matrixDeletion[lineIndex][columnIndex - ONE];
                predecessorValue = valueDeletion;
                break;
            default:
                predecessorCell = matrixInsertion[lineIndex][columnIndex - ONE];
                predecessorValue = valueInsertion;
        }

        matrixInsertion[lineIndex][columnIndex] = new Cell(
                predecessorCell,
                PredecessorType.LEFT,
                firstSequence.charAt(columnIndex - ONE),
                secondSequence.charAt(lineIndex - ONE),
                predecessorValue);
    }

    private void fillDeletionMatrix(int lineIndex, int columnIndex) {
        int valueMatch = matrixMatch[lineIndex - 1][columnIndex].getValue() + OPEN,
                valueInsertion = matrixInsertion[lineIndex - 1][columnIndex].getValue() + OPEN,
                valueDeletion = matrixDeletion[lineIndex - 1][columnIndex].getValue() + EXTEND;
        MatrixType matrixType = MatrixType.DELETION;
        Cell predecessorCell;
        int predecessorValue;

        if (valueMatch >= valueInsertion && valueMatch >= valueDeletion) {
            matrixType = MatrixType.MATCH;
        } else if (valueInsertion >= valueMatch && valueInsertion >= valueDeletion) {
            matrixType = MatrixType.INSERTION;
        }
        switch (matrixType) {
            case MATCH:
                predecessorCell = matrixMatch[lineIndex - ONE][columnIndex];
                predecessorValue = valueMatch;
                break;
            case INSERTION:
                predecessorCell = matrixInsertion[lineIndex - ONE][columnIndex];
                predecessorValue = valueInsertion;
                break;
            default:
                predecessorCell = matrixDeletion[lineIndex - ONE][columnIndex];
                predecessorValue = valueDeletion;
        }

        matrixDeletion[lineIndex][columnIndex] = new Cell(
                predecessorCell,
                PredecessorType.UP,
                firstSequence.charAt(columnIndex - ONE),
                secondSequence.charAt(lineIndex - ONE),
                predecessorValue);
    }

    private void fillScoringMatrix() {
        for (int lineIndex = FIRST_INDEX; lineIndex <= secondSequence.length(); lineIndex++) {
            for (int columnIndex = FIRST_INDEX; columnIndex <= firstSequence.length(); columnIndex++) {
                fillMatchMatrix(lineIndex, columnIndex);
                fillInsertionMatrix(lineIndex, columnIndex);
                fillDeletionMatrix(lineIndex, columnIndex);
            }
        }
    }

    private Cell findMaximumScoreCell() {
        int firstSequenceLength = firstSequence.length(), secondSequenceLength = secondSequence.length();
        int valueMatch = matrixMatch[secondSequenceLength][firstSequenceLength].getValue(),
                valueInsertion = matrixInsertion[secondSequenceLength][firstSequenceLength].getValue(),
                valueDeletion = matrixDeletion[secondSequenceLength][firstSequenceLength].getValue();
        Cell maximumCell;

        if (valueMatch >= valueInsertion && valueMatch >= valueDeletion) {
            maximumCell = matrixMatch[secondSequenceLength][firstSequenceLength];
        } else if (valueInsertion >= valueMatch && valueInsertion >= valueDeletion) {
            maximumCell = matrixInsertion[secondSequenceLength][firstSequenceLength];
        } else maximumCell = matrixDeletion[secondSequenceLength][firstSequenceLength];

        return maximumCell;
    }

    private void align() {
        matrixMatch = new Cell[secondSequence.length() + ONE][firstSequence.length() + ONE];
        matrixInsertion = new Cell[secondSequence.length() + ONE][firstSequence.length() + ONE];
        matrixDeletion = new Cell[secondSequence.length() + ONE][firstSequence.length() + ONE];

        fillGapCells();
        fillScoringMatrix();

        Cell currentCell = findMaximumScoreCell();
        score = currentCell.getValue();
        while (currentCell != null) {
            currentCell = currentCell.fillInformation(firstAlignSequence, secondAlignSequence);
        }
    }

    public PairAlignment(String firstSequence, String secondSequence,
                         ScoringFunction scoringFunction) {
        this.firstSequence = firstSequence;
        this.secondSequence = secondSequence;
        this.scoringFunction = scoringFunction;

        OPEN = getIndel();
        INFINITY = 2 * OPEN + (firstSequence.length() + secondSequence.length()) * EXTEND - 1;

        align();
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        int i = FIRST_INDEX;
        for (; i < Math.ceil((double) firstAlignSequence.length() / NUMBER_OF_SYMBOLS_IN_LINE); i++) {
            stringBuilder.append(OUTPUT_SEQUENCE_1);
            stringBuilder.append(firstAlignSequence.substring(
                    NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE),
                    NUMBER_OF_SYMBOLS_IN_LINE * i));
            stringBuilder.append(NEXT_LINE);

            stringBuilder.append(OUTPUT_SEQUENCE_2);
            stringBuilder.append(secondAlignSequence.substring(
                    NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE),
                    NUMBER_OF_SYMBOLS_IN_LINE * i));
            stringBuilder.append(NEXT_LINE);
            stringBuilder.append(NEXT_LINE);
        }

        stringBuilder.append(OUTPUT_SEQUENCE_1);
        stringBuilder.append(firstAlignSequence.substring(NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE)));
        stringBuilder.append(NEXT_LINE);

        stringBuilder.append(OUTPUT_SEQUENCE_2);
        stringBuilder.append(secondAlignSequence.substring(NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE)));
        stringBuilder.append(NEXT_LINE);
        stringBuilder.append(NEXT_LINE);

        stringBuilder.append(OUTPUT_SCORE);
        stringBuilder.append(score);

        return stringBuilder.toString();
    }
}
