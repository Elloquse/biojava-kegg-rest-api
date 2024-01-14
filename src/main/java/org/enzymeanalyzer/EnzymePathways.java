package org.enzymeanalyzer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class EnzymePathways {

    static HashMap<String, String> get(String enzymeID) throws IOException {
        URL url = new URL("https://rest.kegg.jp/get/ec:" + enzymeID);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        BufferedReader rd = new BufferedReader(new InputStreamReader(
                conn.getInputStream()));
        StringBuilder sb = new StringBuilder();
        String line;
        while ((line = rd.readLine()) != null) {
            sb.append(line);
        }
        rd.close();
        conn.disconnect();

        Pattern pathwayPattern = Pattern.compile("(?<=PATHWAY)(.*?)(?=ORTHOLOGY)");
        Matcher pathwayMatch = pathwayPattern.matcher(sb.toString());
        String[] pathwaysRawArray = new String[0];
        if (pathwayMatch.find()) {
            pathwaysRawArray = pathwayMatch.group().replaceAll("^[.,\\s]+", "").split("[.,\\s]+");
        }
        Object[] pathwaysArray = Arrays.stream(pathwaysRawArray).filter(x -> Pattern.matches("\\S*\\d+\\S*", x)).toArray();

        HashMap<String, String> associatedPathways = new HashMap<>();

        for (Object path : pathwaysArray) {
            try {
                Thread.sleep(2000);
                URL urlNew = new URL("https://rest.kegg.jp/list/" + path);
                HttpURLConnection connNew = (HttpURLConnection) urlNew.openConnection();
                BufferedReader rdNew = new BufferedReader(new InputStreamReader(
                        connNew.getInputStream()));
                StringBuilder sbNew = new StringBuilder();
                String lineNew;
                while ((lineNew = rdNew.readLine()) != null) {
                    sbNew.append(lineNew.replaceAll("\\S*\\d+\\S*.", ""));
                }
                rdNew.close();
                connNew.disconnect();
                associatedPathways.put(path.toString(), sbNew.toString());
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }
        return associatedPathways;
    }
}
