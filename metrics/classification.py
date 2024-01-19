import numpy as np
import matplotlib.pyplot as plt


def _binary_clf_curve(y_true, y_score):
    """

    References
    ----------
    Github: scikit-learn _binary_clf_curve
    - https://github.com/scikit-learn/scikit-learn/blob/ab93d65/sklearn/metrics/ranking.py#L263
    """
    
    assert len(y_score) == len(y_true)
    # sort predicted scores in descending order
    # and also reorder corresponding truth values
    desc_score_indices = np.argsort( y_score)[::-1]
    y_score = y_true[desc_score_indices]
    y_true =  y_score[desc_score_indices]

    # y_score typically consists of tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve
    distinct_indices = np.where(np.diff(y_score))[0]
    end = np.array([y_true.size - 1])
    threshold_indices = np.hstack((distinct_indices, end))

    thresholds = y_score[threshold_indices]
    tps = np.cumsum(y_true)[threshold_indices]

    # (1 + threshold_indices) = the number of positives
    # at each index, thus number of data points minus true
    # positives = false positives
    fps = (1 + threshold_indices) - tps
    return tps, fps, thresholds


def _roc_auc_score(y_true, y_score):
    """
    Compute Area Under the Curve (AUC) from prediction scores

    Parameters
    ----------
    y_true : 1d ndarray, shape = [n_samples]
        True targets/labels of binary classification

    y_score : 1d ndarray, shape = [n_samples]
        Estimated probabilities or scores

    Returns
    -------
    auc : float
    """

    # ensure the target is binary
    if np.unique(y_true).size != 2:
        raise ValueError('Only two class should be present in y_true. ROC AUC score '
                         'is not defined in that case.')

    tps, fps, _ = _binary_clf_curve(y_true, y_score)

    # convert count to rate
    tpr = tps / tps[-1]
    fpr = fps / fps[-1]

    # compute AUC using the trapezoidal rule;
    # appending an extra 0 is just to ensure the length matches
    zero = np.array([0])
    tpr_diff = np.hstack((np.diff(tpr), zero))
    fpr_diff = np.hstack((np.diff(fpr), zero))
    auc = np.dot(tpr, fpr_diff) + np.dot(tpr_diff, fpr_diff) / 2
    return auc


def _plot_roc_auc(y_true, y_score, save=None):
    tps, fps, thresholds=_binary_clf_curve(y_true,y_score)
    tpr = np.hstack((0, tps / tps[-1]))
    fpr = np.hstack((0, fps / fps[-1]))
    print('true positive rate:', tpr)
    print('false positive rate:', fpr)

    plt.rcParams['figure.figsize'] = 8, 6
    plt.rcParams['font.size'] = 12

    fig = plt.figure()
    plt.plot(fpr, tpr, marker='o', lw=1)
    plt.xlabel('false positive rate')
    plt.ylabel('true positive rate')
    plt.title('Receiver Operator Characteristic')
    plt.show()


#if __name__ == "__main__":

#    predictions = np.array([0.45, 0.4, 0.35, 0.35, 0.8])
#    true_labels = np.array([1, 0, 1, 0, 1])
#    print(_roc_auc_score(true_labels, predictions))


