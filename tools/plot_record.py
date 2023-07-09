# > pip install streamlit pandas
# > streamlit run plot_record.py <path/to/log>
#
#
# log is like:
# google::InitGoogleLogging(argv[0]);
# google::InstallFailureSignalHandler();
# google::SetStderrLogging(google::INFO);
# google::SetLogDestination(google::INFO, "tssolver.log.");
#
# LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d})", loop, best_score, score, accept, reject);
import streamlit as st
import pandas as pd
import io
import datetime
import re
import json

pattern = re.compile(r'^(?P<level>.)(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})\s*(?P<hour>\d{2}):(?P<minute>\d{2}):(?P<second>\d{2})\.(?P<microsecond>\d+).*?RECORD (?P<json>.+)$')

uploaded_file = st.file_uploader("Choose a log file")
if uploaded_file is not None:
    fi = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    started_at = None
    rows = []
    for line in fi.readlines():
        mo = pattern.match(line)
        if mo is not None:
            gd = mo.groupdict()
            at = datetime.datetime(
                int(gd['year'], 10),
                int(gd['month'], 10),
                int(gd['day'], 10),
                int(gd['hour'], 10),
                int(gd['minute'], 10),
                int(gd['second'], 10),
                int(gd['microsecond'], 10))
            if started_at is None:
                started_at = at
                td = datetime.timedelta()
            else:
                td = at - started_at

            payload = json.loads(gd['json'])
            payload['time_s'] = td.total_seconds()
            rows.append(payload)
            print(payload)

    df = pd.DataFrame(rows)
    df.set_index('time_s', inplace=True)

    score_columns = []
    if 'best' in df.columns: score_columns.append('best')
    if 'current' in df.columns: score_columns.append('current')
    st.line_chart(df[score_columns])
    if 'accept' in df.columns and 'reject' in df.columns:
        st.line_chart(df[['accept', 'reject']])
    if 'T' in df.columns:
        st.line_chart(df[['T']])

